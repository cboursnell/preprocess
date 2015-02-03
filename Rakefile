require 'rake/testtask'

Rake::TestTask.new do |t|
  t.libs << 'test'
end

Rake::TestTask.new do |t|
  t.name = :trim
  t.libs << 'test'
  t.test_files = ['test/test_test.rb']
end

Rake::TestTask.new do |t|
  t.name = :align
  t.libs << 'test'
  t.test_files = ['test/test_align.rb']
end

Rake::TestTask.new do |t|
  t.name = :normalise
  t.libs << 'test'
  t.test_files = ['test/test_normalise.rb']
end

Rake::TestTask.new do |t|
  t.name = :correction
  t.libs << 'test'
  t.test_files = ['test/test_correction.rb']
end

Rake::TestTask.new do |t|
  t.name = :sra
  t.libs << 'test'
  t.test_files = ['test/test_sra.rb']
end

Rake::TestTask.new do |t|
  t.name = :quant
  t.libs << 'test'
  t.test_files = ['test/test_quantify.rb']
end

desc "Run tests"
task :default => :test