# http://guides.rubygems.org/make-your-own-gem/
class Crimla::Motif
	
	attr_accessor :chr, :start, :stop, :length, :coverage, :consequences, :motif, :type
	
	@@CODONS = {
"TTT" => "Phe",
"TTC" => "Phe",
"TTA" => "Leu",
"TTG" => "Leu",

"TCT" => "Ser",
"TCC" => "Ser",
"TCA" => "Ser",
"TCG" => "Ser",

"TAT" => "Tyr",
"TAC" => "Tyr",
"TAA" => "stop",
"TAG" => "stop",

"TGT" => "Cys",
"TGC" => "Cys",
"TGA" => "stop",
"TGG" => "Trp",

"CTT" => "Leu",
"CTC" => "Leu",
"CTA" => "Leu",
"CTG" => "Leu",

"CCT" => "Pro",
"CCC" => "Pro",
"CCA" => "Pro",
"CCG" => "Pro",

"CAT" => "His",
"CAC" => "His",
"CAA" => "Gln",
"CAG" => "Gln",

"CGT" => "Arg",
"CGC" => "Arg",
"CGA" => "Arg",
"CGG" => "Arg",

"ATT" => "Ile",
"ATC" => "Ile",
"ATA" => "Ile",
"ATG" => "Met",

"ACT" => "Thr",
"ACC" => "Thr",
"ACA" => "Thr",
"ACG" => "Thr",

"AAT" => "Asn",
"AAC" => "Asn",
"AAA" => "Lys",
"AAG" => "Lys",

"AGT" => "Ser",
"AGC" => "Ser",
"AGA" => "Arg",
"AGG" => "Arg",

"GTT" => "Val",
"GTC" => "Val",
"GTA" => "Val",
"GTG" => "Val",

"GCT" => "Ala",
"GCC" => "Ala",
"GCA" => "Ala",
"GCG" => "Ala",

"GAT" => "Asp",
"GAC" => "Asp",
"GAA" => "Glu",

"GAG" => "Gly",
"GGT" => "Gly",
"GGC" => "Gly",
"GGA" => "Gly",
"GGG" => "Gly"
	}
	
	def self.create(opts = {})
		opts = {
			chr: "",
			start: 0,
			stop: 0,
			motif: [],
			exons: [],
			max_inframe: 9,
			refseq: nil
		}.merge(opts)
		# sort the motifs
		opts[:motif].sort!{|m1, m2| m1[:start] <=> m2[:start]}
		motif = Motif.new(opts)
		motif
	end
	
	def initialize(opts)
		opts = {
			chr: "",
			start: 0,
			stop: 0,
			motif: [], 
			exons: [],
			max_inframe: 9,
			refseq: nil
		}.merge(opts)
		@chr = opts[:chr]
		@start = opts[:start]
		@stop = opts[:stop]
		# @start, @stop = [@stop, @start] if @stop < @start
		@length = (opts[:stop] - opts[:start]).abs
		@motif = opts[:motif]
		@consequences = []
		@coverage = [0] * @length 
		@exons = opts[:exons]
		@max_inframe = opts[:max_inframe]
		determine_consequence(opts[:refseq])
	end
	
	def is_wt?()
		@motif.all?{|m| m[:type] == "WT"}
	end
	
	def is_lof?()
		!is_wt? && @consequences.any?{|c| c == "FRAMESHIFT" || c == "INFRAMEDISRUPTION" || c.index("SNV")}
	end
	
	# TODO Calculate consequence of SNV and check if it results in STOP GAIN or STOP LOSS or START GAIN
	# TODO Calculate frameshift using the complete motif. A deletion of 1, followed by a deletion of size 2 might result in an SEMI-INFRAME
	def determine_consequence(refseq = nil)
		if @exons.size > 0 then
			@consequences = @motif.map do |ma|
				consequence = "NA"
				
				if ma[:type] == "WT" then
					consequence = "WILDTYPE"
				else
					mcoords = [ma[:start], ma[:stop]] # this is always 5' -> 3'
					mfrom, mto = mcoords
					
					# determine the length of each exon
					overlaps_exon = @exons.any?{|exstart, exstop| mcoords.any?{|c| c.between?(exstart, exstop) || c.between?(exstop, exstart)}}
					exons_length = @exons.map{|exstart, exstop| (exstop - exstart).abs}
					
					# determine the exon that was hit by the motif and the exons before it
					ehit = @exons.index{|exstart, exstop| mcoords.any?{|c| c.between?(exstart, exstop) || c.between?(exstop, exstart)}}
					exon_hit = @exons[ehit]
					
					reverse_base = {
						"A" => "T",
						"C" => "G",
						"G" => "C",
						"T" => "A",
						"N" => "N"
					}
					
					if ma[:type] == "M" then
						# calculate the reading frame of the SNV
						framehit = nil
						exon_bases = []
						exon_iterator = (0...@exons.length).to_a
						is_reverse = @exons[0][0] > @exons[0][1]
						exon_iterator = exon_iterator.reverse if is_reverse # - strand
						exon_iterator.each do |exoni|
							exfrom, exto = @exons[exoni]
							exonpos = exfrom
							while exonpos != exto do
								exon_base = (refseq[exonpos] || "N")
								exon_base = reverse_base[exon_base] if is_reverse
								exon_bases << exon_base
								if ma[:start] == exonpos then # mark position of SNV
									exon_bases[-1] = "#{exon_bases[-1]}X"
								end
								exonpos += 1 if exfrom < exto
								exonpos -= 1 if exfrom > exto
							end
						end
						
						exon_frames = exon_bases.each_slice(3).to_a
						framehit = exon_frames.index{|frame| frame.any?{|b| b.index("X")}}
						
						if framehit.nil? then
							consequence = "SNV(#{ma[:ref]}>#{ma[:seq]})"
						else
							# replace reference with mutated base and determine amino acid
							mutated_pos = exon_frames[framehit].index{|b| b.index("X")}
							codon_seq = exon_frames[framehit].dup
							if !is_reverse then
								codon_seq[mutated_pos] = ma[:seq]
							else
								codon_seq[mutated_pos] = reverse_base[ma[:seq]]
							end
							codon_seq = codon_seq.join("").gsub("X", "")
							
							codon_ref = exon_frames[framehit].join("").gsub("X", "")
							
							aref = @@CODONS[codon_ref].to_s
							aseq = @@CODONS[codon_seq].to_s
							if aref == "" then # in this case we don't know the refernce of the coodrinates were not given correctly
								STDERR.puts "[WARNING] Could not determine reference codon at #{ma[:start]}"
								consequence = "NOREF(#{ma[:ref]}>#{ma[:seq]})"
							elsif aref == aseq then
								consequence = "SYN(#{ma[:ref]}>#{ma[:seq]})"
							else
								if aseq == "stop" then
									consequence = "STOPGAIN(#{ma[:ref]}>#{ma[:seq]})"
								elsif aref == "stop"
									consequence = "STOPLOSS(#{ma[:ref]}>#{ma[:seq]})"
								elsif aseq == "Met"
									consequence = "STARTGAIN(#{ma[:ref]}>#{ma[:seq]})"
								elsif aref == "Met"
									consequence = "STARTLOSS(#{ma[:ref]}>#{ma[:seq]})"
								else
									consequence = "SNV(#{ma[:ref]}>#{ma[:seq]}|#{aref}>#{aseq})"
								end
							end
							consequence
						end
						
					else
						
						if !overlaps_exon then # INTRON, we have to do nothing.
							consequence = "INTRON"
						else
							
							# calculate the number of nucleotides before the motif
							# this depends on the direction of the translation
							if @exons[0][0] < @exons[0][1] then # + strand
								exon_length_before_hit = exons_length[0...ehit]
								nucleotides_before_hit = (exon_length_before_hit + [(mfrom - exon_hit[0]).abs]).inject{|sum,x| sum + x }
							else # - strand
								exon_length_before_hit = exons_length[(ehit+1)..-1]
								nucleotides_before_hit = (exon_length_before_hit + [(mto - exon_hit[0]).abs]).inject{|sum,x| sum + x }
							end
							
							motif_length = mto - mfrom
							
							# check for INFRAME or FRAMESHIFT mutations
							if nucleotides_before_hit % 3 == 0 and motif_length % 3 == 0 then # we have a hit inside the reading frame
								if motif_length < @max_inframe || @max_inframe < 0
									consequence = "INFRAME"
								else
									consequence = "INFRAMEDISRUPTION"
								end
							else
								consequence = "FRAMESHIFT"
							end
						end
					end
				end # end of if type == WT
				consequence
			end # end of map
		end # end of if @exons.size > 0
	end
	
	def increase_coverage(from, to, by = 1, offset = @start)
		deletions = @motif.select{|ma| ma[:type] == "D"}.map{|ma| (ma[:start]...ma[:stop])}
		((from-offset)...(to-offset)).each do |i|
			next if i >= @coverage.size or i < 0
			next if deletions.any?{|d| d.cover?(i+offset)}
			@coverage[i] += by
		end
	end
	
	def to_s()
		# "#{@chr}:#{@start}-#{@stop}/" + @motif.each_with_index.map{|m, i| "#{m[:type]}:#{m[:start]}-#{m[:stop]}[#{@consequences[i]}]"}.join("/")
		# "#{@chr}:#{@start}-#{@stop}/" + @motif.each_with_index.map{|m, i| "#{m[:start]}:#{m[:type]}:#{(m[:start]-m[:stop]).abs}[#{@consequences[i]}]"}.join("/")
		@motif.each_with_index.map{|m, i| "#{m[:type]}_#{(m[:start]-m[:stop]).abs}bp:#{m[:start]}-#{m[:stop]}[#{@consequences[i]}]"}.join("/")
	end
end