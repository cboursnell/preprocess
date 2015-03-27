#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestNormalise < Test::Unit::TestCase

  context 'normalise' do

    setup do
      input = File.join(File.dirname(__FILE__), 'data', 'raw_data')
      @output = Dir.mktmpdir
      verbose = false
      threads = 1
      memory = 1
      @pre = Preprocessor::Preprocessor.new(@output, verbose, threads, memory)
      @pre.load_input(input)
    end

    teardown do
      # delete output folder
      cmd = "rm -rf #{@output}"
      `#{cmd}`
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      bindir = File.join(gem_dir, "bin")
      cmd = "rm #{bindir}/trimmomatic-0.32.jar"
      `#{cmd}` if File.exist?(File.join(bindir, "trimmomatic-0.32.jar"))
      cmd = "rm -rf #{bindir}/bbmap"
      `#{cmd}` if Dir.exist?(File.join(bindir, "bbmap"))
      cmd = "rm #{bindir}/bbmap.tar.gz"
      `#{cmd}` if File.exist?(File.join(bindir, "bbmap.tar.gz"))
    end

    should 'run khmer' do
      @pre.trimmomatic
      @pre.khmer(23, 5, 4)
      assert File.exist?("#{@output}/test.khmered.fq-left.fq")
      assert File.exist?("#{@output}/test.khmered.fq-right.fq")
      assert File.exist?("#{@output}/test.unpaired.fq")
    end

    should 'run bbnorm' do
      @pre.bbnorm
      assert File.exist?("#{@output}/test-A-1_1.bbnorm.fq"), "file exist 1-1"
      assert File.exist?("#{@output}/test-A-1_2.bbnorm.fq"), "file exist 1-2"
      assert File.exist?("#{@output}/test-A-2_1.bbnorm.fq"), "file exist 2-1"
      assert File.exist?("#{@output}/test-A-2_2.bbnorm.fq"), "file exist 2-2"
    end

    should 'trim and then bbnorm' do
      @pre.trimmomatic
      @pre.bbnorm
      assert File.exist?("#{@output}/test-A-1_1.bbnorm.fq"), "file exist 1-1"
      assert File.exist?("#{@output}/test-A-1_2.bbnorm.fq"), "file exist 1-2"
      assert File.exist?("#{@output}/test-A-2_1.bbnorm.fq"), "file exist 2-1"
      assert File.exist?("#{@output}/test-A-2_2.bbnorm.fq"), "file exist 2-2"
    end

  end
end
