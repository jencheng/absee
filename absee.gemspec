Gem::Specification.new do |s|
  s.name        = 'absee'
  s.version     = '1.1.1'
  s.date        = '2018-12-19'
  s.summary     = ".ab1 reader / ABIF reader"
  s.description = ".ab1 reader / ABIF reader; extracts the peak indexes, called sequence, quality scores, and ACGT values from sequencing files"
  s.authors     = ["Jenny Cheng"]
  s.email       = 'jencheng@ginkgobioworks.com'
  s.files       = ["lib/absee.rb"]
  s.required_ruby_version = '>= 1.9.3'
  s.license     = 'MIT'
  s.homepage    = 'http://rubygems.org/gems/absee'
  s.metadata = {
    "documentation_uri" => "https://github.com/jencheng/absee"
  }
end
