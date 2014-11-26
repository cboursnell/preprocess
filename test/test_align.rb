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

    should 'run bowtie2 with single reads' do
      reference = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      pre = Preprocessor::Preprocessor.new(@output, false, 1, 1)
      left = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq')
      left << ","
      left << File.join(File.dirname(__FILE__), 'data', 'A-2-1.fq')
      pre.load_single_reads(left, "test")
      pre.bowtie2(reference, false)
      assert File.exist?("#{@output}/A-1-1-reference.fa.sam")
      assert File.exist?("#{@output}/A-2-1-reference.fa.sam")
      assert File.exist?("#{@output}/bowtie2.stats")
    end

    should 'run snap' do
      reference = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      @pre.snap(reference)
      assert File.exist?("#{@output}/A-1-1.fq.A-1-2.fq.reference.bam")
      assert File.exist?("#{@output}/A-2-1.fq.A-2-2.fq.reference.bam")
      assert File.exist?("#{@output}/snap.stats")
    end

  end
end
