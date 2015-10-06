# hold many motifs
require "motif"
class MotifPopulation
	
	attr_accessor :chr, :start, :stop, :length, :coverage, :consequences, :motif, :type, :exons
	
	def initialize(opts)
		opts = {
			chr: "",
			start: 0,
			stop: 0,
			exons: []
		}.merge(opts)
		@MOTIFS = {}
		@MOTIF_COUNT = Hash.new(0)
		@chr = opts[:chr]
		@start = opts[:start]
		@stop = opts[:stop]
		# @start, @stop = [@stop, @start] if @stop < @start
		@length = (opts[:stop] - opts[:start]).abs
		@coverage = [0] * @length 
		@exons = opts[:exons]
	end
	
	def size
		@MOTIF_COUNT.size
	end
	
	def self.reset()
		@MOTIFS = {}
		@MOTIF_COUNT = Hash.new(0)
	end
	
	def apply_count_threshold(threshold)
		@MOTIF_COUNT.delete_if{|motif, cnt| cnt < threshold}
		@MOTIF_COUNT
	end
	
	def find_or_create(opts = {})
		opts = {
			chr: "",
			start: 0,
			stop: 0,
			motif: [],
			exons: [],
			max_inframe: 9
		}.merge(opts)
		# sort the motifs
		opts[:motif].sort!{|m1, m2| m1[:start] <=> m2[:start]}
		motif = Motif.new(opts)
		if @MOTIFS[opts].nil? then
			@MOTIFS[opts] = motif
		else
			motif = @MOTIFS[opts]
		end
		motif
	end
	
	def increase_count(motif, by = 1)
		@MOTIF_COUNT[motif] += by
		@MOTIF_COUNT[motif]
	end
	
	def top_motifs(top = nil, types = %w(D I M WT), &block)
		sorted_motifs = @MOTIF_COUNT.keys.sort{|m1, m2| @MOTIF_COUNT[m2] <=> @MOTIF_COUNT[m1]}.select{|m| m.motif.any?{|ma| types.include?(ma[:type])} }
		total_reads = sum(@MOTIF_COUNT.values.flatten)
		if !top.nil? then
			top = sorted_motifs.size if top > sorted_motifs.size
			sorted_motifs = sorted_motifs[0...top]
		end
		total_top_reads = sum(sorted_motifs.map{|m| @MOTIF_COUNT[m]})
		
		curr_cumfreq = 0
		curr_top_cumfreq = 0
		curr_count_lof = 0
		curr_rank = 0
		sorted_motifs.map!{|m|
			total_freq = @MOTIF_COUNT[m].to_f / total_reads.to_f
			total_top_freq = @MOTIF_COUNT[m].to_f / total_top_reads.to_f
			curr_cumfreq += total_freq
			curr_top_cumfreq += total_top_freq
			curr_rank += 1
			if m.is_lof?() then
				curr_count_lof += @MOTIF_COUNT[m].to_f
			end
			curr_lof_cumfreq = curr_count_lof / total_reads.to_f
			curr_top_lof_cumfreq = curr_count_lof / total_top_reads.to_f
			{
				motif: m, 
				rank: curr_rank,
				count: @MOTIF_COUNT[m], 
				freq_total: total_freq,
				freq_total_top: @MOTIF_COUNT[m].to_f / total_top_reads,
				cumfreq: curr_cumfreq,
				cumfreq_top: curr_top_cumfreq,
				cumfreq_lof: curr_lof_cumfreq,
				cumfreq_top_lof: curr_top_lof_cumfreq
			}
		}
		
		if block_given? then
			sorted_motifs.each do |m|
				yield m
			end
		else
			return sorted_motifs
		end
	end
	
	def purity_index(top = 2)
		return -1 if top < 1
		top_count = sum(top_motifs(top, %w(D I M WT)).map{|m| m[:count]})
		total_support = sum(@MOTIF_COUNT.values.flatten)
		top_count.to_f / total_support.to_f
	end
	
	# find out how many motifs are required to explain a proportion of reads
	def find_quantiles(q = [50, 90, 95, 99])
		topm = top_motifs(nil, %w(D I M WT))
		ret = {}
		q.each{|qv|
			ret[qv] = (topm.find_index{|m|
				m[:cumfreq] >= qv.to_f/100
			}.to_f + 1).to_f
		}
		ret
	end
	
	def entropy(top = nil, base = 2)
		rel_motif = top_motifs(top).map{|m| m[:freq_total_top]}
		
		log_method = :log2
		log_method = :log if base == :e or base == "e"
		log_method = :log10 if base == 10
		h = -1 * sum(
			rel_motif.map{|p| p * Math.send(log_method, p)}
		)
		return h
	end
	
	def trapezoidArea(x1, y1, x2, y2)
		return ((y1 + y2)/2.0) * (x2 - x1)
	end
	
	# field should be cumfreq, cumfreq_top, cumfreq_lof or cumfreq_top_lof
	def auc(top = nil, field = :cumfreq_top_lof)
		points = top_motifs(top).map{|m| [m[:rank], m[field]]}
		# now scale the rank
		points.map!{|p| [p[0].to_f/points.size.to_f, p[1]]}
		# add artificial starting point
		points = [[0,0]] + points
		
		calculate_auc(points)
	end
	
	def calculate_auc(points)
		traps = (1..(points.size-1)).map{|i|
			p = points[i]
			px = points[i-1]
			ret = trapezoidArea(px[0], px[1], p[0], p[1])
		}
		auc = sum(traps)
		auc
	end
	
	def num_motifs(types = %w(D I WT), top = nil)
		mytypes = Hash[types.map{|t| [t, true]}]
		mytypes.default = false
		if top.nil? then
			@MOTIF_COUNT.keys.select{|m| m.motif.any?{|ma| mytypes[ma[:type]]} }.size
		else
			mymotifs = self.top_motifs(top).map{|m| m[:motif]}
			mymotifs.select{|m| m.motif.any?{|ma| mytypes[ma[:type]]} }.size
		end
	end
	
	
	def create_histograms(top = nil, types = %w(D I))
		motifs = self.top_motifs(top, types)
		hists = {}
		types.each{|t| 
			hists[t] = {
				start: Hash.new(0),
				stop: Hash.new(0),
				length: Hash.new(0),
				region: Hash.new(0)
			}
		}
		motifs.each do |motif_count|
			motif = motif_count[:motif]
			count = motif_count[:count]
			motif.motif.each do |ma|
				next if hists[ma[:type]].nil?
				hists[ma[:type]][:start][ma[:start]] += count
				hists[ma[:type]][:stop][ma[:stop]] += count
				hists[ma[:type]][:length][(ma[:start] - ma[:stop]).abs] += count
				(ma[:start]...ma[:stop]).each do |i|
					hists[ma[:type]][:region][i] += count
				end
			end
		end
		hists
	end
	
	def cumulate_coverage(start, stop, top = @MOTIFS.size)
		cumulated_coverage = []
		return [] if @MOTIFS.size == 0
		chromosome = @MOTIFS[@MOTIFS.keys.first].chr
		self.top_motifs(top).each do |motif_cnt|
			motif = motif_cnt[:motif]
			raise "all motifs have to be on the same chromosome" if motif.chr != chromosome
			motif.coverage.each_with_index do |coverage, i|
				next if i+motif.start-start < 0 
				next if i+motif.start-start > (stop-start)
				cumulated_coverage[i+motif.start-start] = 0 if cumulated_coverage[i+motif.start-start].nil?
				cumulated_coverage[i+motif.start-start] += coverage
			end
		end
		cumulated_coverage
	end
	
private
	def sum(a)
		a.inject{|sum,x| sum + x }.to_f
	end
end