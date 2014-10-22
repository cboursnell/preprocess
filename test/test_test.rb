#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestPreprocessor < Test::Unit::TestCase

  context 'preprocessor' do

    setup do
      input = File.join(File.dirname(__FILE__), 'data', 'raw_data')
      @output = Dir.mktmpdir
      @verbose = false
      @threads = 1
      @memory = 1
      @pre = Preprocessor::Preprocessor.new(@output, @verbose, @threads, @memory)
      @pre.load_input(input)
    end

    teardown do
      # delete output folder
      cmd = "rm -rf #{@output}"
      `#{cmd}`
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      bindir = File.join(gem_dir, "bin")
      trimmo_path = File.join(ENV['GEM_HOME'], 'bin', 'trimmomatic-0.32.jar')
      # cmd = "rm #{trimmo_path}"
      # `#{cmd}` if File.exist?(trimmo_path)
      cmd = "rm -rf #{bindir}/SPAdes-3.1.0-Linux"
      `#{cmd}` if Dir.exist?(File.join(bindir, "SPAdes-3.1.0-Linux"))
      cmd = "rm #{bindir}/spades.tar.gz"
      `#{cmd}` if File.exist?(File.join(bindir, "spades.tar.gz"))
    end

    should 'setup should run ok' do
      assert @pre
    end

    should 'install trimmomatic' do
      @pre.trimmomatic
      paths = which('trimmomatic-0.32.jar')
      if paths.empty?
        path = File.join(ENV['GEM_HOME'], 'bin', 'trimmomatic-0.32.jar')
      else
        path = paths.first
      end
      cmd = "java -jar #{path}"
      java_cmd = Preprocessor::Cmd.new(cmd)
      java_cmd.run
      assert_equal 357, java_cmd.stderr.length
      assert_equal "Usage: ", java_cmd.stderr.split("\n").first
    end

    should 'trim reads using trimmomatic' do
      @pre.trimmomatic
      assert File.exist?("#{@output}/A_1-1.t.fq"), "file doesn't exist"
      assert File.exist?("#{@output}/A_1-2.t.fq"), "file doesn't exist"
      @pre.data.each do |hash|
        assert hash[:current]
        assert hash[:unpaired]
      end
    end

    should 'trim single reads using trimmomatic' do
      @pre = Preprocessor::Preprocessor.new(@output, @verbose, @threads, @memory)
      reads = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq')
      reads << ",#{File.join(File.dirname(__FILE__), 'data', 'A-2-1.fq')}"
      @pre.load_single_reads(reads, "A")
      @pre.trimmomatic
      assert File.exist?("#{@output}/A_1.t.fq"), "A_1.t.fq doesn't exist"
      assert File.exist?("#{@output}/A_2.t.fq"), "A_2.q.fq doesn't exist"
    end

    # should 'install skewer' do
    #   @pre.skewer
    #   cmd = "skewer-0.1.117-linux-x86_64 --help"
    #   skewer_cmd = Preprocessor::Cmd.new(cmd)
    #   skewer_cmd.run
    #   stdout = skewer_cmd.stdout.split("\n")
    #   str="Skewer (A fast and accurate adapter trimmer for paired-end reads)"
    #   ver="Version 0.1.117 (updated in July 12, 2014), Author: Hongshan Jiang"
    #   assert_equal str, stdout[0]
    #   assert_equal ver, stdout[1]
    # end

    should 'trim reads using skewer' do
      @pre.skewer
      @pre.data.each do |hash|
        assert File.exist?(hash[:current]), "file exist"
        assert hash[:current]=~/test-A-[12]-pair/, "file name"
      end
    end

    should 'load reads without input file' do
      verbose = false
      left = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq')
      left << ",#{File.join(File.dirname(__FILE__), 'data', 'A-2-1.fq')}"
      right = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq')
      right << ",#{File.join(File.dirname(__FILE__), 'data', 'A-2-2.fq')}"
      pre = Preprocessor::Preprocessor.new(@output, verbose, @threads, @memory)
      pre.load_reads(left, right, "A")
      pre.trimmomatic
      pre.hammer
      pre.bbnorm
    end

    should "gunzip files and leave original in place" do
      left = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq.gz')
      right = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq.gz')
      output = Dir.mktmpdir
      pre = Preprocessor::Preprocessor.new(output, false, 1, 1)
      pre.load_reads(left, right, "A")
      pre.gunzip
      gz1 = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq.gz')
      gz2 = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq.gz')
      assert File.exist?(gz1), "original file 1 exists"
      assert File.exist?(gz2), "original file 2 exists"
      assert File.exist?("#{output}/A-1-1.fq"), "gunzipped file 1 exists"
      assert File.exist?("#{output}/A-1-2.fq"), "gunzipped file 2 exists"
      assert_equal "#{output}/A-1-1.fq", pre.data[0][:current]
      assert_equal "#{output}/A-1-2.fq", pre.data[1][:current]
    end

  end
end
