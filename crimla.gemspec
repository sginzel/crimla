Gem::Specification.new do |s|
	s.name        = 'crimla'
	s.version     = '0.0.1'
	s.date        = '2015-08-20'
	s.summary     = "A tool to analyse the diversity of genome modifications."
	s.description = "This tool can be used to analyse the diversity of CRISPR experiments."
	s.authors     = ["Sebastian Ginzel"]
	s.email       = 'sginze2s@inf.h-brs.de'
	s.files       = ["lib/crimla.rb", "lib/crimla/motif.rb", "lib/crimla/motif_population.rb", "lib/crimla/report.rb", "lib/crimla/parser.rb", "bin/crimla", "LICENSE"]
	s.homepage    = 'http://github.com/sginzel/crimla'
	s.license     = 'MIT'
	s.executables << 'crimla'
end