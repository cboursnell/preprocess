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
      pre.load_reads(left, right)
      pre.trimmomatic
      pre.hammer
      pre.bbnorm
    end

    should 'trim reads using trimmomatic' do
      @pre.trimmomatic
      assert File.exist?("#{@output}/A_1-1.t.fq")
      assert File.exist?("#{@output}/A_1-2.t.fq")
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
