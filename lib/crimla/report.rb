class Crimla::Report
	
	def initialize(motif_population, parser)
		@population = motif_population
		@parser = parser
	end
	
	def calculate_stats(opts)
		count = @parser.count
		
		uniq_motifs = @population.num_motifs
		uniq_motif_per_pair = uniq_motifs.to_f / count[:pairs]
		total_uniq_motif_with_lof = @population.top_motifs(nil, %w(D I)).select{|ma| ma[:motif].is_lof?}.size
		top_uniq_motif_with_lof = @population.top_motifs(opts[:top], %w(D I)).select{|ma| ma[:motif].is_lof?}.size
		fraction_pairs_indel = @population.top_motifs(nil, %w(D I)).select{|ma| !ma[:motif].is_wt?}.map{|ma| ma[:count]}.inject{|sum,x| sum + x }.to_f / count[:pairs].to_f
		fraction_pairs_lof = @population.top_motifs(nil, %w(D I)).select{|ma| ma[:motif].is_lof?}.map{|ma| ma[:count]}.inject{|sum,x| sum + x }.to_f / count[:pairs].to_f
		
		auc_diversity = @population.auc(nil, :cumfreq)
		auc_lof = @population.auc(nil, :cumfreq_lof)
		entropy = @population.entropy
		shannon = @population.entropy(nil, "e")
		eveness = shannon / Math.log(uniq_motifs).to_f
		diversity = Math.exp(shannon)
		quantiles = @population.find_quantiles([25, 75, 90, 95, 99])
		iqr = quantiles[25]/quantiles[75]
		spread = 1-iqr.to_f
		spread = (quantiles[99]-quantiles[90]).to_f/quantiles[99]
		q90 = quantiles[90]/uniq_motifs
		q95 = quantiles[95]/uniq_motifs
		q99 = quantiles[99]/uniq_motifs
		clonal_purity_index = @population.purity_index(opts["purity-ploidy"])
		
		@stats = {
			uniq_motifs: uniq_motifs,
			uniq_motif_per_pair: uniq_motif_per_pair,
			total_uniq_motif_with_lof: total_uniq_motif_with_lof,
			top_uniq_motif_with_lof: top_uniq_motif_with_lof,
			fraction_pairs_indel: fraction_pairs_indel,
			fraction_pairs_lof: fraction_pairs_lof,
			auc_diversity: auc_diversity,
			auc_lof: auc_lof,
			entropy: entropy,
			diversity: diversity,
			shannon: shannon,
			eveness: eveness,
			quantiles: quantiles,
			iqr: iqr,
			spread: spread,
			clonal_purity_index: clonal_purity_index,
			q90: q90,
			q95: q95,
			q99: q99
		}
		
	end
	
	def add_text_block(name, p, start, stop, labels = [], fill_color = "white")
		block_track = p.add_track(
			:glyph => :label, 
			:name => name, 
			:label => true, 
			:fill_color => fill_color,
			:track_height=>10
		)
		block = Bio::Graphics::MiniFeature.new(
			:start => (start-(stop-start)*0.1).round(0), 
			:end => (stop+(stop-start)*0.1).round(0), 
			:id => nil,
			label: false
		)
		block_track.add(block)
		
		labels.each do |str|
			offset = str.scan(/^[ ]*/).first.size * 2
			line = Bio::Graphics::MiniFeature.new(
				:start => start + offset, 
				:end => start + offset + 1, 
				:id => str.to_s,
				label: false
			)
			block_track.add(line)
		end
		block_track
	end
	
	def generate(opts)
		
		calculate_stats(opts)
		
		chr = @population.chr
		start = @population.start
		stop = @population.stop
		exons = @population.exons
		count = @parser.count
		
		if !opts[:svg].nil? then
		################# SVG ########################
		STDERR.puts "Generating SVG..."
		p = Bio::Graphics::Page.new(
			:width => 1200,
			:height => 200,
			:number_of_intervals => 8
		)
		add_text_block("Statistics", p, start, stop, labels = [
			"Cell Population",
			"  Distinct motifs: #{@population.num_motifs}",
			"  Distinct motifs/K#{(opts[:paired])?"pair":"read"}: #{(@stats[:uniq_motif_per_pair]*1000).round(4)}",
			"  Diversity: #{@stats[:diversity].round(3)}",
#			"  Spread: #{@stats[:spread].round(3)}", 
			"  Purity index: #{@stats[:clonal_purity_index].round(3)}",
			"  Shannon: #{@stats[:shannon].round(3)}",
			"  Eveness: #{@stats[:eveness].round(3)}",
			"KOE", 
			"  Unique motifs with LOF: #{@stats[:total_uniq_motif_with_lof]}",
			"  #{(opts[:paired])?"Pairs":"Reads"} supporting INDEL-motifs: #{@stats[:fraction_pairs_indel].round(3)*100}%",
			"  #{(opts[:paired])?"Pairs":"Reads"} with LOF: #{@stats[:fraction_pairs_lof].round(3)*100}%",
			"Statistical Measures",
			"  Entropy: #{@stats[:entropy].round(3)}",
			"  Quantiles: #{@stats[:quantiles].map{|q, r| "Q#{q}=#{r}"}.join(", ")}",
			"  Rel. quantiles: Q90=#{@stats[:q90].round(3)}, Q95=#{@stats[:q95].round(3)}, Q99=#{@stats[:q99].round(3)}"
			
		], fill_color = "white")
		
		if exons.size > 0 then
			gene_track = p.add_track(
				:glyph => :transcript, 
				:name => "Exons", 
				:label => true, 
				:exon_fill_color => "grey",
				:exon_stroke_width => 0,
				:exon_style => "fill-opacity:0.5;",
				:track_height => 200,
				:height => 30
			)
			exon_feature =	Bio::Graphics::MiniFeature.new(
					:start => start, 
					:end => stop, 
					:exons => exons.flatten,
					:id => nil
				)
			gene_track.add(exon_feature)
		end # end of add exons
		
		if opts[:coverage] && opts["coverage-global"]
			global_coverage = @population.cumulate_coverage(start, stop, nil)#opts[:top])
			track_gobal_coverage = p.add_track(
				:glyph => :histogram,
				:track_height => 100,
				:stroke_width => 0,
				:start => (start-(stop-start)*0.1).round(0), 
				:end => (stop+(stop-start)*0.1).round(0), 
				:fill_color => "grey",
				:label => true, 
				:name => "Global Region Coverage of #{@population.num_motifs} motifs (max: #{global_coverage.max})"
			)
			global_coverage.each_with_index do |coverage, i|
				track_gobal_coverage.add(Bio::Graphics::MiniFeature.new(
					:start => (start+i),
					:end => (start+i+1),
					:segment_height => coverage
				))
			end
		end # of gloabal coverage
		
		# we add one to the top parameter because we show reference also
		@population.top_motifs((opts[:top].nil?)?nil:(opts[:top].to_i), %w(D I WT)).each do |motif_count|
			motif_counter = motif_count[:rank]
			motif = motif_count[:motif]
			mcount = motif_count[:count].to_f
			cum_freq = motif_count[:cumfreq].to_f.round(3) * 100
			freq = motif_count[:freq_total].round(3)*100
			
			motifstr = motif.to_s
			if motif.motif.any?{|ma| ma[:type] == "WT"} then
				motif_type = "Wildtype"
				fill_color = "green"
				motifprop = (mcount/count[:wt].to_f).round(3)*100
				delexons = motif.motif.select{|ma| ma[:type] == "WT"}.map{|ma| [ma[:start], ma[:stop]]}.flatten
				insexons = []
				styles = Hash.new("fill-opacity:1.0;")
			else
				motif_type = "Motif"
				fill_color = "blue"
				styles = {0 => "fill-opacity:0.6;", 1 => "fill-opacity:0.3;"}
				motifprop = (mcount/count[:motif].to_f).round(3)*100
				# insexons = motif.motif.select{|ma| ma[:type] == "I"}.map{|ma| [ma[:start], ma[:stop]]}.flatten
				insexons = motif.motif.select{|ma| ma[:type] == "I"}.map{|ma| [ma[:start], ma[:start]+0.5]}
				snpexons = motif.motif.select{|ma| ma[:type] == "M"}.map{|ma| [ma[:start], ma[:start]+1]}
				tmp = []
				tmp += insexons unless insexons.size == 0
				tmp += snpexons unless snpexons.size == 0
				# tmp += snpexons.reverse if insexons.size == 0 # this is wierd, but this way SNPs are encoded as block and inserions as arrows.
				insexons = tmp.sort{|x,y| x[0] <=> y[0]}.flatten(1) # snps and insersion are both represented the same way, but should be sorted
				
				delexons = motif.motif.select{|ma| ma[:type] == "D"}.map{|ma| [ma[:start], ma[:stop]]}.flatten
				if delexons.size > 0 then
					delexons = [start] + delexons + [stop]
				else
					delexons = [start, stop]
				end
				insexons = [start, start] + [stop, stop] + insexons ## this hack removes wierd blocks from the ends of the exons...
			end
			#cum_freq += freq
			motif_track = p.add_track(
				:glyph => :transcript, 
				:name => "##{motif_counter}) #{motif_type}: #{motifstr}, Freq: #{freq.round(3)}%, Cum.Freq: #{(cum_freq).round(3)}%, #{motif_type}.Freq: #{motifprop.round(3)}%, Abs. count: #{mcount.to_i}",#"#{m.pretty_inspect}", 
				:label => true, 
				:exon_fill_color => fill_color, #"blue",
				:exon_stroke_width => 0,
				:exon_style => styles[motif_counter % 2], # "fill-opacity:0.4;",
				:utr_fill_color => (insexons.size > 4)?"red":"blue", 
				:utr_stroke_width => 0,
				:utr_style => (insexons.size > 4)?"fill-opacity:1.0;":"fill-opacity:0.1;",
				:x_round => 0,
				:y_round => 0
			)
			motif_feature =	Bio::Graphics::MiniFeature.new(
					:start => start, 
					:end => stop, 
					:exons => delexons,
					:utrs => insexons,
					:id => nil
				)
			motif_track.add(motif_feature)
			# Add Histogram
			if opts[:coverage] && opts["coverage-motif"] then
				hist_track = p.add_track(
					:glyph => :histogram,
					:track_height => 50,
					:stroke_width => 0,
					:label => true, 
					:fill_color => fill_color,
					:style => styles[motif_counter % 2],
					:name => "Coverage (max: #{motif.coverage.max})"
				)
				motif.coverage.each_with_index do |coverage, i|
					hist_track.add(Bio::Graphics::MiniFeature.new(
						:start => (start+i),
						:end => (start+i+1),
						:segment_height => coverage
					))
				end
			end
		end
		
		## add histogram of deletion/insertion distribution
		hists = @population.create_histograms(nil, %w(D I))
		[
			["D", {name: "Deletion (start/stop) [green/orange]", fields: [:start, :stop], color: "orange"}],
			["I", {name: "Insertion (start/stop) [green/orange]", fields: [:start, :stop], color: "orange"}],
			["D", {name: "Deletion region", fields: [:region], color: "black"}],
			["I", {name: "Insertion region", fields: [:region], color: "black"}]
		].each do |motif_type, conf|
			max_start = (max = hists[motif_type][:start].values.flatten.max; hists[motif_type][:start].keys.select{|pos| hists[motif_type][:start][pos] == max}.join(", "))
			max_stop  = (max = hists[motif_type][:stop].values.flatten.max; hists[motif_type][:stop].keys.select{|pos| hists[motif_type][:stop][pos] == max}.join(", "))
			max_region  = (max = hists[motif_type][:region].values.flatten.max; hists[motif_type][:region].keys.select{|pos| hists[motif_type][:region][pos] == max}.join(", "))
			maxes = {
				start: max_start,
				stop: max_stop,
				region: max_region
			}
			hist_track = p.add_track(
					:glyph => :histogram,
					:track_height => 50,
					:stroke_width => 0,
					:label => true, 
					:fill_color => conf[:color],
					:style => "fill-opacity:0.5;",
					:name => "#{conf[:name]} (freq. #{conf[:fields].map{|f| "#{f} pos.: " + maxes[f].to_s}.join(", ")})"
				)
			conf[:fields].each do |field|
				hists[motif_type][field].each do |pos, poscount|
					hist_track.add(Bio::Graphics::MiniFeature.new(
					{
						:start => (pos),
						:end => (pos+1),
						:params => {
							#fill_color: (field == :start)?"dark#{conf[:color]}":"#{conf[:color]}"
							fill_color: (field == :start)?"green":"#{conf[:color]}"
						},
						:segment_height => poscount
					}))
				end
			end
		end
		
		add_text_block("Processing Statistics", p, start, stop, labels = [
			"Reads", 
			"  Total (Reads; Pairs): #{count[:total]}; #{count[:pairs]}, Processed: #{count[:parsed]}",
			"  Supporting WT/Motifs: #{count[:wt]}/#{count[:motif]}",
			"Filters", 
			"  Low.Qual: #{count[:skipped_qual]}, Too short: #{count[:skipped_len]}", 
			"  Off-Region #{count[:off_target]}, Skipped (other): #{count[:skipped]}", 
			"  No-Overlap: #{count["no-overlap"]}, Not-processed: #{count["not-processed"]}",
			"  Low. Coverage: #{count[:low_coverage]}"
		], fill_color = "white")
		
		## add parameters to SVG
		param_track = p.add_track(
			:glyph => :label, 
			:name => 'Parameter', 
			:label => true, 
			:track_height=>10
		)
		opts.each do |k,v|
			next if k.to_s =~ /_given$/
			postfix = ""
			if opts["#{k}_given".to_sym].nil? then
				postfix = " (default)"
			end
			param_text = Bio::Graphics::MiniFeature.new(
				:start => start, 
				:end => start + 1, 
				:id => "#{k}: #{v}#{postfix}",
				label: false
			)
			param_track.add(param_text)
		end
		
		## add LEGEND
		legend_track = p.add_track(
			:glyph => :label, 
			:name => 'Legend', 
			:label => true, 
			:track_height=>10
		)
		
		{
			"SECTION: Statistics" => "",
			"Cell Population" => " ",
			"Distinct motifs" => "Number of distinct motifs (including Wildtype WT)",
			"Distinct motifs/K(read|pair)" => "Ratio of distinct motifs per 1000 reads or read-pairs (if paired). Range [0,1]",
			"Diversity" => "Diversity Index defined as exp(-sum(pm * ln pm)), where pm is the proporion of reads/read-pairs supporting motif m. [0,] - higher value equals higher diversity.",
#			"Spread" => "Def: ((Q99-Q90)/Q99), Description: Ratio of motifs that explain 90% of the data. Use Q99 instead of Q100 to remove outliers. Range [0,1] - low spread to high spread",
			"Purity index" => "Ratio of support for the two most abundant motifs and the reads/read-pairs for all motifs. The two most abundant motifs should be equal to the two most abundant alleles, if the clone is pure. Range [0, 1]: low purity - high purity",
			"Shannon index" => "Shannon index of population defined as (-sum(pm * ln pm), where pm is the proporion of reads/read-pairs supporting motif m. [0,] - higher value equals higher diversity.",
			"Eveness" => "Shannon index divided by maximum possible entropy. [0,1] - higher value equals higher eveness.",
			"KOE (knock-out-efficiency)" => " ",
			"Unique motifs with LOF" => "Number of distinct motifs with frameshift INDEL",
			"Reads/Pairs supporting INDEL-motifs" => "Ratio of motifs containing at least one INDEL.",
			"Reads/Pairs with LOF" => "Ratio of motif containing at least one INDEL that leads to predicted LOF",
			"Statistical Measures" => " ",
			"Entropy" => "Shannon index for motif counts. Defined as -sum(pm * log2 pm)), where pm is the proporion of reads/read-pairs supporting motif m. [0,] - lower values mean less abundance",
			"Quantiles" => "Number of motifs required to explain 25%, 75%, 90%, 95%, or 99% of the observed reads/read-pairs. Range [0,] - the number of motifs is not scaled or bound.",
			"Rel. quantiles" => "Numer of motifs required to explain a set proportion of the data, divided by the total number of motifs. Range [0,1] - higher value means more data is explained.",
			"SECTION: Global Region Coverage" => "",
			"Description" => "Accumulated coverage of all reads. Overlapping part of read pairs contribute 2X coverage.",
			"SECTION: Motifs" => "",
			"Motif" => "Description of the observed motif, which can include many submotifs. First letters indicate type of submotif (I D WT), first number is the length followd by the start and stop coordinates.",
			"Freq" => "Frequency of the motif over all reads",
			"Cum.Freq" => "Cumulated frequency of this motif and all motifs which are more abundant.",
			"Motif/Wildtype Freq." => "Fraction of the occurence of the motif compared to all other motifs of the same class (WT or INDEL)",
			"Abs. count" => "Number of pairs supporting this motif",
			"SECTION Deletion/Insertion/Region" => "",
			"start/stop" => "Histogram of the location of indel start and stop locations for all submotifs.",
			"freq. start pos / max. stop pos" => "Location of the most frequent start and stop positions(peaks in the histogram)",
			"Region" => "Histogram of the regions affected by deletions or insertions including the location of the peak.",
			"SECTION: Processing Statistics" => "", 
			"Reads" => " ",
			"Total (Reads; Pairs)" => "Total number of reads/read-pairs in the BAM file",
			"Processed" => "Number of reads/read-pairs parsed (after filtering)",
			"WT / Motif" => "Number of reads/read-pairs supporting wildtype/not wildtype allele",
			"Filters" => " ",
			"Low.Qual / Too short" => "Number of reads which do not match the quality or length threshold. Hard clipped parts do not contribute to read length.",
			"Off-Region / Skipped" => "Number of reads that are not inside the query region / Number of reads skipped due to inconsistencies (such as both reads of a pair being tagged as the first in pair)",
			"No overlap / Not-processed" => "Number of read pairs that do not overlap or are singletons",
			"SECTION: Parameter" => "",
			"List of parameters used to create the SVG" => "Mostly used for reproduction of the image.",
			"SECTION: Legend" => "",
			"Legend" => "This legend."
		}.each do |k,v|
			offset = (v.length == 0)?0:7
			offset -= (v == " ")?3:0
			legend_text = Bio::Graphics::MiniFeature.new(
				:start => start + offset, 
				:end => start + offset + 1, 
				:id => "#{k}#{(v.length == 0)?"":":"} #{v}",
				label: false
			)
			legend_track.add(legend_text)
		end
		
		STDERR.print "writing to #{opts[:svg]}..."
		begin
			p.write(opts[:svg])
			STDERR.print " DONE\n"
		rescue => e
			STDERR.print e.message
			STDERR.print " WARNING COULD NOT WRITE SVG - continuing to write stats...\n"
		end
		
		
	end # of svg generation
	
	if !opts[:csv].nil? then
		STDERR.puts "writing stats to #{opts[:csv]}"
		fout = File.new(opts[:csv], "w+")
	else
		fout = STDOUT
	end
	
	recs = []
	# cum_freq = 0
	num_frameshifts = 0
	rank = 0
	@population.top_motifs(nil, %w(D I WT)).each_with_index do |motif_count, i|
		rank = i + 1
		num_frameshifts += (motif_count[:motif].is_lof?)?1:0
		# freq = motif_count[:count].to_f/count[:pairs].to_f
		# cum_freq += freq
		recs << motif_count.merge({
			is_wt: motif_count[:motif].is_wt?,
			num_lof: num_frameshifts,
			total_pairs: count[:pairs],
			total_wt: count[:wt],
			total_motif: count[:motif],
			start_motif: motif_count[:motif].motif.map{|ma| ma[:start]}.min,
			end_motif: motif_count[:motif].motif.map{|ma| ma[:stop]}.max
		}).merge(@stats)
		#recs << {
		#	rank: rank,
		#	motif: motif_count[:motif].to_s, 
		#	count: motif_count[:count],
		#	is_wt: motif_count[:motif].is_wt?,
		#	freq: (freq).round(3),
		#	cumulated_freq: cum_freq.round(3),
		#	total_pairs: count[:pairs],
		#	total_wt: count[:wt],
		#	total_motif: count[:motif],
		#	num_lof: num_frameshifts,
		#	freq_lof: (num_frameshifts.to_f/rank.to_f).round(3)
		#}
	end
	fout.puts (recs.first || {}).keys.join("\t")
	puts (recs.first || {}).keys.join("\t") unless opts[:csv].nil?
	recs.each do |rec|
		fout.puts rec.values.join("\t")
		puts rec.values.join("\t") unless opts[:csv].nil?
	end
	end
	
end