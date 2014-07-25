#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestPreprocessor < Test::Unit::TestCase

  context 'preprocessor' do

    setup do
      input = File.join(File.dirname(__FILE__), 'data', 'raw_data')
      @output = Dir.mktmpdir
      verbose = false
      threads = 1
      memory = 1
      @pre = Preprocessor.new(@output, verbose, threads, memory)
      @pre.load_input(input)
    end

    teardown do
      # delete output folder
      cmd = "rm -rf #{@output}"
      `#{cmd}`
    end

    should 'setup should run ok' do
      assert @pre
    end

    should 'load reads without input file' do
      verbose = false
      threads = 1
      memory = 1
      left = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq')
      left << ",#{File.join(File.dirname(__FILE__), 'data', 'A-2-1.fq')}"
      right = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq')
      right << ",#{File.join(File.dirname(__FILE__), 'data', 'A-2-2.fq')}"
      pre = Preprocessor.new(@output, verbose, threads, memory)
      pre.load_reads(left, right, "A")
      pre.trimmomatic
      pre.hammer
      pre.bbnorm
    end

    should "gunzip files and leave original in place" do
      left = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq.gz')
      right = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq.gz')
      output = Dir.mktmpdir
      pre = Preprocessor.new(output, false, 1, 1)
      pre.load_reads(left, right, "A")
      pre.gunzip
      gz1 = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq.gz')
      gz2 = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq.gz')
      assert File.exist?(gz1), "original file 1 exists"
      assert File.exist?(gz2), "original file 2 exists"
      assert File.exist?("#{output}/A-1-1.fq"), "gunzipped file 1 exists"
      assert File.exist?("#{output}/A-1-2.fq"), "gunzipped file 2 exists"
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

    should 'normalise reads with khmer' do
      @pre.khmer
      assert File.exist?("#{@output}/A.left.fq")
      assert File.exist?("#{@output}/A.right.fq")
    end

    should 'construct hammer input' do
      @pre.construct_hammer_input
      assert File.exist?("#{@output}/dataset.yaml")
    end

    should 'hammer reads' do
      @pre.hammer
    end

    should 'run bbnorm' do
      @pre.bbnorm
      assert File.exist?("#{@output}/A-1-1.bbnorm.fq")
      assert File.exist?("#{@output}/A-1-2.bbnorm.fq")
      assert File.exist?("#{@output}/A-2-1.bbnorm.fq")
      assert File.exist?("#{@output}/A-2-2.bbnorm.fq")
    end

    should 'trim and then bbnorm' do
      @pre.trimmomatic
      @pre.bbnorm
      assert File.exist?("#{@output}/A-1-1U.bbnorm.fq")
      assert File.exist?("#{@output}/A-1-2U.bbnorm.fq")
      assert File.exist?("#{@output}/A-2-1U.bbnorm.fq")
      assert File.exist?("#{@output}/A-2-2U.bbnorm.fq")
    end

    # should 'run my normaliser' do
    #   @pre.norm
    # end
  end
end
