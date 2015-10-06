Gem::Specification.new do |s|
	s.name        = 'ramplicon'
	s.version     = '0.0.1'
	s.date        = '2015-08-20'
	s.summary     = "A tool to analyse the diversity of genome modifications."
	s.description = "This tool can be used to analyse the diversity of CRISPR experiments."
	s.authors     = ["Sebastian Ginzel"]
	s.email       = 'sginze2s@inf.h-brs.de'
	s.files       = ["lib/motif.rb", "lib/motif_population.rb", "lib/report.rb", "lib/parser.rb", "ramplicon"]
	s.homepage    = 'http://rubygems.org/gems/ramplicon'
	s.license     = 'MIT'
	s.executables << 'ramplicon'
end