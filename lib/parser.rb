require "report"
require "motif"
require "motif_population"

# Class used to 

class Parser
	
	@@cigarex = Regexp.new(/(([0-9]+)([MSID]))/) # split cigar by operands - this does not include N and P operands. H operand is also not included as it does not contribute to the coordinates.
	
	attr_accessor :count
	
	def initialize()
		@count = {total: 0, pairs: 0, skipped_qual: 0, skipped_len: 0, skipped: 0, parsed: 0, "no-overlap" => 0, wt: 0, motif: 0, off_target: 0}
		@count.default = 0
	end
	
	def increase_count(field)
		@count[field] += 1
	end
	
	# returns a motif population
	def process(opts)
		STDERR.puts "#"*80
		STDERR.puts "# Processing #{opts[:bam]}..."
		STDERR.puts "#"*80
		opts.each do |k,v|
			STDERR.puts "   #{k}: #{v}"
		end
		STDERR.puts "#"*80
		
		chr, start, stop = opts[:coord].to_s.split(/[:-]/).map{|x| x.gsub(",", "")}
		start = start.to_i
		stop = stop.to_i
		if stop < start then
			start, stop = [stop, start]
		end
		
		start = start - opts[:tail].to_i
		stop  = stop + opts[:tail].to_i
		
		exons = opts[:exons]
							.map{|e| 
								elems = e.gsub(",", "").split(/[-:]/, 3)
								if elems.size == 3 then
									if elems[0] != chr then
										STDERR.puts "Exon and Coordinates are not on the same chromosome"
										exit 1
									end
									elems[1..2].map(&:to_i) # first one is chromosome
								else
									elems[0..1].map(&:to_i)
								end
								}.uniq
							.sort{|e1, e2| e1[0] <=> e2[0]}
		
		# adjust the length of the exon that includes to start of the reading frame
		if !opts[:atg].nil? then
			atgchr, atg = opts[:atg].gsub(",", "").split(":", 2)
			atg = atgchr if atg.nil?
			atg = atg.to_i
			exons.map!{|estart, estop|
				if atg.between?(estart, estop) || atg.between?(estop, estart) then
					estart = atg
				end
				[estart, estop]
			}
		end
		
		mp = MotifPopulation.new(chr: chr, start: start, stop: stop, exons: exons)
		
		cnt = 0
		@count = {total: 0, pairs: 0, skipped_qual: 0, skipped_len: 0, skipped: 0, parsed: 0, "no-overlap" => 0, wt: 0, motif: 0, off_target: 0}
		@count.default = 0
		
		pairs = Hash.new(nil)
		incomplete_pairs = 0
		
		num_paired_reads = 1 # reads required to start looking for motifs...
		num_paired_reads = 2 if opts[:paired] 
		
		if opts[:qnames].to_a.size > 0 then
			qnames_given = opts[:qnames].map{|qn|
				if File.exists?(qn)
					File.new(qn, "r").readlines.select{|l| l[0] != "#"}.map(&:strip)
				else 
					qn
				end
			}.flatten.uniq
			qnames = Hash.new(0)
			qnames_given.each do |qn|
				qnames[qn] = (opts[:paired])?2:1
			end
			qnames_cnt = qnames.size
		else
			qnames = Hash.new(0)
			qnames_cnt = 0
		end
		
		STDERR.puts "opening BAM for coordinates #{chr}:#{start}-#{stop}"
		bam = Bio::DB::Sam.new(:bam => opts[:bam], :fasta => opts[:ref])
		if !bam.indexed? then
			STDERR.puts "Indexing BAM file..."
			bam.index
		end
		
		# MyDeBug = Struct.new(:pos, :cigar, :mapq, :seq)
		motifs = Hash.new(0)
		qnames_to_skip = Hash.new(false)
		start_at = Time.now
		reads_per_second = 0
		reads_per_minute = 0
		
		refseq = nil
		if opts[:snvs] then
			refseq = bam.fetch_reference("#{chr}", start.to_i, stop.to_i).to_s 
			if refseq == "" then
				STDERR.puts "Reference sequence not found."
				exit 2
			end
			# make reference a Hash with the corresponding coordinates
			refseq = Hash[refseq.split("").each_with_index.map{|base, i| [start+i, base]}]
		end
		
		bam.fetch("#{chr}", start.to_i, stop.to_i) do |a|
		# [MyDeBug.new(27950262, "48M20D51M", 100, "N"*119)].each do |a|
			if cnt % 10000 == 0
				reads_per_second = (cnt.to_f / ((Time.now - start_at).to_f).to_f)
				reads_per_minute = (reads_per_second*60).round(0).to_i
				reads_per_minute = reads_per_minute.to_s.reverse.gsub(/(\d{3})(?=\d)/, '\\1,').reverse
				STDERR.print "[#{cnt}/#{@count[:parsed]}/#{@count[:pairs]}/#{pairs.size}]#{(opts[:skip] > cnt)?"skipping":"processing"} (#{reads_per_minute} rpm)...                       \r"
			end 
			cnt += 1
			
			next if opts[:skip] > cnt
			
			qname = a.qname
			cigar = a.cigar
			strand = a.query_strand ? '+' : '-'
			pairno = a.first_in_pair ? 0 : 1
			# pairno = (strand == "+") ? 0 : 1
			paired_size = a.isize
			
			# if there are qnames then check fo valid qname
			if opts[:qnames].size > 0 then
				break if qnames_cnt == 0 # there are no more read names to process...
				next unless qnames[a.qname] > 0
				qnames_cnt -= 1 if (qnames[a.qname] -= 1) == 0
			end
			
			increase_count(:total)
			
			if !a.pos.between?(start, stop) || !(a.pos+a.seq.length).between?(start, stop) then
				increase_count(:off_target)
				qnames_to_skip[a.qname] = true
				pairs.delete(a.qname) #unless pairs[a.qname].nil? # delete the pair if the partner has a bad quality
				next
			end
			
			if (opts["stop-after"].to_i > 0) then
				if @count[:pairs] >= opts["stop-after"].to_i
					STDERR.puts "Stopping after #{@count[:pairs]} processed pairs..."
					break
				end 
			end
			
			if (opts[:paired]) then
				if not (a.is_paired && a.is_mapped && !a.mate_unmapped) then
					increase_count(:skipped)
					next 
				end
			end
			
			if a.mapq.to_i < opts["mapping-quality"].to_i
				increase_count(:skipped_qual)
				qnames_to_skip[a.qname] = true
				pairs.delete(qname) # delete the pair if the partner has a bad quality
				next
			end 
			
			# perform filtering and statistics for each read
			# but skip further processing if it was already black listed.
			# we can't do this earlier because it get inconsistent with paired reads.
			# when proccessing single end reads we will not skip them because reads should be
			# proccessed and judged individually
			next if qnames_to_skip[qname] && opts[:paired]
			
			## make an entry for a new pair
			if pairs[qname].nil? then
				pairs[qname] = {
					paired_size: paired_size,
					nreads: 0,
					reads: []
				}
			end
			# add current read to its pair
			if !pairs[qname][:reads][pairno].nil? then
				# in this case the first_in_pair attribute was not set correctly and we skip the read pair
				increase_count(:skipped)
				qnames_to_skip[a.qname] = true
				pairs.delete(qname)
				next
			end
			
			if opts[:snvs] then
				aseq = a.seq.split("").each_with_index.map{|b,i| [b, a.qual[i].ord - 33]}
			else
				aseq = nil
			end 
			
			pairs[qname][:reads][pairno] = {
				pos: a.pos,
				cigar: cigar,
				strand: strand,
				seq: aseq#(opts[:snvs])?(a.seq):("")
			}
			pairs[qname][:nreads] += 1
			
			# if there are enough reads for the read group start to look into the reads for motifs
			if pairs[qname][:nreads] == num_paired_reads then
				## FIND motifs
				# puts "found pair #{qname}"
				
				# check for inconsistencies. This can happen when both reads are marked as "first in pair"
				pairs[qname][:reads].reject!(&:nil?)
				if pairs[qname][:reads].size != pairs[qname][:nreads] then
					increase_count(:skipped)
					qnames_to_skip[a.qname] = true # make sure that we skip reads with the same name - although there shouldn't be more paired reads with the same ID
					pairs.delete(qname)
					next
				end
				
				min_support = num_paired_reads
				if !opts["overlap-only"] then
					min_support = 0
				end
				
				pairs[qname][:motifs] = find_motifs(pairs[qname][:reads], min_support, refseq, opts["min-base-quality"].to_i, opts["min-overlap"].to_i, start, stop)
				
				increase_count(:parsed)
				
				# After determining the motif we can determine the length of the reads
				if pairs[qname][:reads].any?{|r| r[:length] < opts["min-len"].to_i} then
					increase_count(:skipped_len)
					qnames_to_skip[a.qname] = true
					pairs.delete(qname)
					next
				end
				
				# find_motifs return nil if the reads didnt overlap to determine the motif
				# it returns an empty motif if too few reads support a motif
				if pairs[qname][:motifs].nil? or pairs[qname][:motifs].size == 0 then 
					# @count["no-overlap"] += 1
					increase_count("no-overlap")
					qnames_to_skip[a.qname] = true 
					pairs.delete(qname)
					next
				else
					motifs[pairs[qname][:motifs]] += 1 # this was used to @count before the Motif.increase_@count method
				end
				
				# id the pair is successful and overlaps 
				increase_count(:pairs)
				# get the motif as a singleton
				motif = mp.find_or_create({
					chr: chr,
					start: start,
					stop: stop,
					motif: pairs[qname][:motifs].keys,
					exons: exons,
					refseq: refseq,
					rfoffset: opts["reading-frame-offset"],
					max_inframe: opts["tolerated-indel-len"]
				})
				mp.increase_count(motif, 1)
				# @count motifs
				if motif.is_wt? then
					# @count[:wt] += 1
					increase_count(:wt)
				else
					# @count[:motif] += 1
					increase_count(:motif)
				end
				# add coverage from supporting reads to the motif...
				if opts[:coverage] || opts["min-motif-count"] > 0 then
					pairs[qname][:reads].map{|r| [r[:pos], r[:endpos]]}.sort{|r1,r2| r1[0] <=> r2[0]}.each do |r|
						motif.increase_coverage(r[0], r[1])
					end
				end
				pairs.delete(qname) # we are done with the reads and can discard the partners.
			end
		end # of bam.fetch
		STDERR.print "[#{cnt}/#{@count[:parsed]}/#{@count[:pairs]}/#{pairs.size}]#{(opts[:skip] > cnt)?"skipping":"processing"} (#{reads_per_minute} rpm)...                       \r"
		STDERR.print "DONE parsing the bam file\n"
		@count["not-processed"] = pairs.size
		@count[:low_coverage] = mp.size
		STDERR.print "Filtering #{@count[:low_coverage]} motifs with coverage lower than #{opts["min-motif-count"]}..."
		mp.apply_count_threshold(opts["min-motif-count"]) if opts["min-motif-count"] > 0
		@count[:low_coverage] = @count[:low_coverage] - mp.size
		STDERR.print "DONE.\n"
		mp
	end
	
	# returns a Hash of motifs and how many of the given reads support them
	# if min_support > 0 then only motifs that are supported by more than min_support reads are returned 
	# {
	# 	{ :start=>27950330, 
	#			:stop=>27950334, 
	#			:len=>4, 
	#			:type=>"D"} => 2,
	#		{ :start=>27950340, 
	#			:stop=>27950344, 
	#			:len=>4, 
	#			:type=>"I"} => 2, 
	#		{ :start=>27950267, 
	#			:stop=>27950268, 
	#			:len=>1, 
	#			:type=>"D"} => 1
	#	}
	def find_motifs(reads, min_support = 0, refseq = nil, min_base_qual = 30, min_overlap = -10, from, to)
		
		motifs = Hash.new(0) # create a hash that contains the motif as well as the number of times it was covered
		reads.each do |r|
			next if r.nil?
			cops = r[:cigar].scan(@@cigarex) # extract cigar operators
			currpos = r[:pos]
			refpos = r[:pos]
			seqpos = 0
			readlen = 0
			cops.each do |op, len, type|
				# process SNVs and INDELS
				if type == "M" and !refseq.nil? then
					read_range = (refpos...(refpos+len.to_i)) #if r[:strand] == "+"
					#read_range = ((currpos+len.to_i)..currpos) if r[:strand] == "-"
					read_range.each_with_index{|refi, seqi|
						refb = refseq[refi]
						seqb = r[:seq][seqi+seqpos][0]
						seqb_qual = r[:seq][seqi+seqpos][1]
						if refb != seqb and seqb_qual >= min_base_qual then
							motif = {
								start: refi,
								stop: refi,
								len: 0,
								type: "M",
								ref: refb,
								seq: seqb 
							}
							motifs[motif] += 1
						end
					}
				elsif type == "D" or type == "I"
					motif = {
						start: refpos,
						stop: refpos + len.to_i,
						len: len.to_i,
						type: type
					}
					motifs[motif] += 1
				end
				
				if type == "M" or type == "D" 
					#currpos.send(direction, len.to_i)
					refpos += len.to_i
				end
				if type == "M" or type == "I"
					readlen += len.to_i
					seqpos += len.to_i
				end
			end
			r[:endpos] = refpos - 1#((r[:strand] == "-")?1:0)
			r[:length] = readlen
		end
		
		# check if the reads overlap
		reads_overlap = true
		if reads.reject(&:nil?).size > 1 then # if the read size is 1 there is no overlap to check
			raise "Cannot parse multi-read files" if reads.reject(&:nil?).size > 2
			sorted_reads = reads.map{|r| [r[:pos], r[:endpos]]}.sort{|r1,r2| r1[0] <=> r2[0]} 
			num_overlap_bases = (sorted_reads[1][0] - sorted_reads[0][1]) * -1
			reads_overlap = num_overlap_bases >= min_overlap
			# reads_overlap = reads[0][:pos].between?(reads[1][:pos], reads[1][:endpos]) || reads[1][:pos].between?(reads[0][:pos], reads[0][:endpos])
		end
		return nil if not reads_overlap # if the reads dont overlap we will discard the motif they produce
		
		# if there are no motifs we assume WT - we already checked if the reads overlap so we can assume WT for the region
		if motifs.size == 0 then
			# reads.each{|r| r[:pos] -= 25} # TODO remove this
			start = reads.map{|r| r[:pos]}.min
			stop  = reads.map{|r| r[:endpos]}.max
			motif = {
				# start: start, 
				# stop: stop,
				# len: (stop+1) - start,
				start: from,
				stop: to,
				len: (to+1) - from,
				type: "WT"
			}
			motifs[motif] = reads.size
		end
		# remove motifs that are not supported by both reads
		if min_support > 0 then
			motifs.select!{|k,v| v >= min_support}
		end
		motifs
	end
	
end