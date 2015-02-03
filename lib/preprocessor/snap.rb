require 'bindeps'
require 'csv'

module Preprocessor

  class Snap

    def initialize(outdir, threads, reference, expression)
      @outdir = outdir
      @threads = threads
      @reference = File.expand_path(reference)
      @expression = expression

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'snap.yaml')
      Bindeps.require gem_deps
      @snap = get_bin_path("snap")
      @stats={:reads=>0,:paired_reads=>0,:unaligned=>0,:unique=>0,:multi=>0}
    end

    def run(left, right)
      build_index(@reference, @threads)
      map_reads(left, right, @threads)
    end

    def get_bin_path bin
      which_bin = Cmd.new("which #{bin}")
      which_bin.run
      if !which_bin.status.success?
        raise IOError.new("Snap: could not find #{bin} in path")
      end
      which_bin.stdout.split("\n").first
    end

    def build_paired_cmd left, right, threads
      cmd = "#{@snap} paired #{@outdir}/#{@index_name}"
      cmd << " #{left[:current]}"
      cmd << " #{right[:current]}"
      cmd << " -o #{@outdir}/#{@bam}"
      cmd << " -s 0 1000" # min and max distance between paired-read starts
      cmd << " -H 300000" # max seed hits to consider in paired mode
      cmd << " -h 2000" # max seed hits to consider when reverting to single
      cmd << " -d 30" # max edit distance (function of read length?)
      cmd << " -t #{threads}"
      cmd << " -b" # bind threads to cores
      cmd << " -M"  # format cigar string
      if @expression
        cmd << " -om 5" # Output multiple alignments. extra edit distance
        cmd << " -omax 20" # max alignments per pair/read
      end
      cmd
    end

    def build_single_cmd reads, threads
      cmd = "#{@snap} single #{@outdir}/#{@index_name}"
      cmd << " #{reads[:current]}"
      cmd << " -o #{@outdir}/#{@bam}"
      cmd << " -H 300000" # max seed hits to consider in paired mode
      cmd << " -h 2000" # max seed hits to consider when reverting to single
      cmd << " -d 30" # max edit distance (function of read length?)
      cmd << " -t #{threads}"
      cmd << " -b" # bind threads to cores
      cmd << " -M"  # format cigar string
      if @expression
        cmd << " -om 5" # Output multiple alignments. extra edit distance
        cmd << " -omax 20" # max alignments per pair/read
      end
      cmd
    end

    def map_reads(left, right, threads)
      raise RuntimeError.new("Snap index not built") if !@index_built

      lbase = File.basename(left[:current])
      rbase = File.basename(right[:current])
      @bam = "#{lbase}.#{rbase}.#{@index_name}.bam"

      unless File.exists? "#{@outdir}/#{@bam}"
        snapcmd = build_paired_cmd(left, right, threads)
        runner = Cmd.new snapcmd
        runner.run
        @stats = runner.stdout
        unless runner.status.success?
          raise RuntimeError.new("Snap failed\n#{runner.stderr}")
        end
      end
      left[:sam] = @bam
      right[:sam] = @bam
    end

    def build_index file, threads
      @index_name = File.basename(file, File.extname(file))
      unless Dir.exists?("#{@outdir}/#{@index_name}")
        cmd = "#{@snap} index #{file} #{@outdir}/#{@index_name}"
        cmd << " -s 20"
        cmd << " -t#{threads}"
        cmd << " -bSpace" # contig name terminates with space char
        runner = Cmd.new cmd
        runner.run
        if !runner.status.success?
          err = runner.stderr
          msg = "Failed to build Snap index\n#{runner.stderr}"
          raise RuntimeError.new(msg)
        end
      end
      @index_built = true
    end

    def get_stats
      return @stats
    end

  end

end
