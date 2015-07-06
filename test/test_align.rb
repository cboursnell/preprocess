#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestAlign < Test::Unit::TestCase

  context 'aligner' do

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
    end

    should 'run bowtie2' do
      reference = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      @pre.trimmomatic
      @pre.bowtie2(reference, false)
      assert File.exist?("#{@output}/A_1-1.t-A_1-2.t-reference.fa.sam")
      assert File.exist?("#{@output}/A_2-1.t-A_2-2.t-reference.fa.sam")
      assert File.exist?("#{@output}/bowtie2.stats")
    end

    should 'run bwa index' do
      ref = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      bwa = Preprocessor::Bwa.new(@output, 1, ref, false)
      index = bwa.build_index
      assert File.exist?("#{index}.bwt"), "index exists"
    end

    should 'run bwa mem' do
      ref = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      left = {}
      right = {}
      left[:current] = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq')
      right[:current] = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq')
      bwa = Preprocessor::Bwa.new(@output, 1, ref, false)
      bwa.run left, right
      assert File.exist?(left[:sam])
    end

    # TODO:
    # should 'run snap with single reads' do
    # end

  end
end
