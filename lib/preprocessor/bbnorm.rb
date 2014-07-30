require 'which'
require 'bindeps'
include Which

module Preprocessor

  class BBnorm

    def initialize(outdir, threads, memory, k, coverage, bits, tables,
                   lowthresh, mindepth, minkmers)
      @outdir = outdir
      @threads = threads
      @memory = memory
      @k = k
      @coverage = coverage
      @bits = bits
      @tables = tables
      @lowthresh = lowthresh
      @mindepth = mindepth
      @minkmers = minkmers

      @java = which('java').first

      @gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      @bbnorm = File.join(@gem_dir, 'bin', 'bbmap', 'bbnorm.sh')
      if !File.exist?(@bbnorm)
        install
      end
    end

    def install
      dl = "#{@gem_dir}/bin/bbmap.tar.gz"
      cmd = "wget http://kent.dl.sourceforge.net/project/bbmap/"
      cmd << "BBMap_33.08_java7.tar.gz -O #{dl}"
      wget = Cmd.new(cmd)
      puts wget.cmd
      wget.run
      puts wget.stdout
      if wget.status.success?
        cmd = "tar xzf #{dl} --directory "
        cmd << "#{@gem_dir}/bin"
        untar = Cmd.new(cmd)
        untar.run
        puts untar.stdout
        puts untar.cmd
        if !untar.status.success?
          raise RuntimeError.new("untar failed")
        end
      else
        raise RuntimeError.new("download failed")
      end
    end

    def run left, right=nil
      if left and right
        cmd = "#{@bbnorm} "
        cmd << "in=#{left[:current]} "
        cmd << "in2=#{right[:current]} "
        name = "#{left[:name]}-#{left[:type]}-#{left[:rep]}.bbnorm.fq"
        outfile = File.join(@outdir, name) # interleaved output of left and right
        cmd << "out=#{outfile} "
        cmd << "k=#{@k} "
        cmd << "hashes=#{@tables} "
        cmd << "bits=#{@bits} "
        cmd << "lowthresh=#{@lowthresh} "
        cmd << "mindepth=#{@mindepth} "
        cmd << "minkmers=#{@minkmers} "
        cmd << "target=#{@coverage} "
        cmd << "threads=#{@threads}"
        norm = Cmd.new(cmd)
        norm.run
        if !norm.status.success?
          msg = "something went wrong with bbnorm\n"
          msg << "#{norm.stdout}\n#{norm.stderr}"
          raise RuntimeError.new(msg)
        end
        leftoutput = "#{@outdir}/"
        leftoutput << "#{left[:name]}-#{left[:type]}-"
        leftoutput << "#{left[:rep]}-#{left[:pair]}.bbnorm.fq"
        rightoutput = "#{@outdir}/"
        rightoutput << "#{left[:name]}-#{right[:type]}-"
        rightoutput << "#{right[:rep]}-#{right[:pair]}.bbnorm.fq"
        self.deinterleave(outfile, leftoutput, rightoutput)
        left[:current] = leftoutput
        right[:current] = rightoutput
      end
    end

    def deinterleave(file, output_left, output_right)
      #puts "deinterleaving #{file} to make #{output_left} and #{output_right}"
      raise RuntimeError.new("Can't find #{file}") if !File.exist?(file)
      fastq = File.open(file)
      left = File.open("#{output_left}", "w")
      right= File.open("#{output_right}", "w")

      name1 = fastq.readline rescue nil
      seq1 = fastq.readline rescue nil
      plus1 = fastq.readline rescue nil
      qual1 = fastq.readline rescue nil
      name2 = fastq.readline rescue nil
      seq2 = fastq.readline rescue nil
      plus2 = fastq.readline rescue nil
      qual2 = fastq.readline rescue nil

      while name1 != nil and name2 != nil
        left.write(name1)
        left.write(seq1)
        left.write(plus1)
        left.write(qual1)
        right.write(name2)
        right.write(seq2)
        right.write(plus2)
        right.write(qual2)

        name1 = fastq.readline rescue nil
        seq1 = fastq.readline rescue nil
        plus1 = fastq.readline rescue nil
        qual1 = fastq.readline rescue nil
        name2 = fastq.readline rescue nil
        seq2 = fastq.readline rescue nil
        plus2 = fastq.readline rescue nil
        qual2 = fastq.readline rescue nil
      end
      fastq.close
      left.close
      right.close
    end

  end

end