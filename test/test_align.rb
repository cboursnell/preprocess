#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestAlign < Test::Unit::TestCase

  context 'aligner' do

    setup do
      input = File.join(File.dirname(__FILE__), 'data', 'raw_data')
      @output = Dir.mktmpdir
      puts @output
      verbose = true
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
    end

    should 'run bowtie2' do
      reference = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      @pre.trimmomatic
      @pre.bowtie2(reference, false)
      assert File.exist?("#{@output}/A_1-1.t-A_1-2.t-reference.fa.sam")
      assert File.exist?("#{@output}/A_2-1.t-A_2-2.t-reference.fa.sam")
      assert File.exist?("#{@output}/bowtie2.stats")
    end

    should 'run bowtie2 and express' do
      reference = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      @pre.trimmomatic
      @pre.bowtie2(reference, true)
      @pre.express(reference)
      assert File.exist?(File.join(@output, "express_results"))
      assert File.exist?(File.join(@output, "test_DE_ebseq.txt"))


    end

  end
end
