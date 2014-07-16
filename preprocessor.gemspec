Gem::Specification.new do |gem|
  gem.name          = 'preprocessor'
  gem.version       = '0.1'
  gem.date          = '2014-05-27'
  gem.summary       = "Preprocess mRNA reads"
  gem.description   = "See summary"
  gem.authors       = ["Chris Boursnell", "Richard Smith-Unna"]
  gem.email         = 'cmb211@cam.ac.uk'
  gem.files         = `git ls-files`.split("\n")
  gem.executables   = ["preprocess"]
  gem.require_paths = %w( lib )
  gem.homepage      = 'http://rubygems.org/gems/preprocessor'
  gem.license       = 'MIT'

  gem.add_dependency 'which', '~> 0.0', '>= 0.0.2'
  gem.add_dependency 'bio', '~> 1.4', '>= 1.4.3'
  gem.add_dependency 'bindeps', '~> 0.0', '>= 0.0.7'

  gem.add_development_dependency 'rake', '~> 10.3', '>= 10.3.2'
  gem.add_development_dependency 'turn', '~> 0.9', '>= 0.9.7'
  gem.add_development_dependency 'simplecov', '~> 0.8', '>= 0.8.2'
  gem.add_development_dependency 'shoulda-context', '~> 1.2', '>= 1.2.1'
  gem.add_development_dependency 'coveralls', '~> 0.7'
end