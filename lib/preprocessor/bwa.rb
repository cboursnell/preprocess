require 'bindeps'
require 'csv'
require 'fixwhich'

module Preprocessor

  class Bwa

    def initialize(outdir, threads, reference, expression)
      @outdir = outdir
      @threads = threads
      @reference = File.expand_path(reference)
      @expression = expression

      # gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      # gem_deps = File.join(gem_dir, 'deps', 'bwa.yaml')
      # Bindeps.require gem_deps
      @bwa = Which.which('bwa').first
      # @stats={:reads=>0,:paired_reads=>0,:unaligned=>0,:unique=>0,:multi=>0}
    end

    # Usage:   bwa index [options] <in.fasta>
    def build_index
      # construct index if it doesn't exist
      index = File.join(@outdir, "#{File.basename(@reference)}")
      if !File.exist?("#{index}.bwt")
        cmd = "#{@bwa} index"
        cmd << " #{@reference}"
        build = Cmd.new(cmd)
        build.run
        if !build.status.success?
          msg = "Something went wrong with bwa index\n"
          msg << "#{build.stdout}"
          msg << "#{build.stderr}"
          raise RuntimeError.new(msg)
        end
        ["bwt", "pac", "ann", "amb", "sa"].each do |ext|
          File.rename("#{@reference}.#{ext}", "#{index}.#{ext}")
        end
      end
      return index
    end

    # Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
    def run left, right=nil
      index = build_index
      l = File.basename(left[:current])
      r = File.basename(right[:current]) if right
      sam = "#{l.split(".")[0..-2].join(".")}-"
      sam << "#{r.split(".")[0..-2].join(".")}-" if right
      sam << "#{File.basename(index)}.sam"
      sam = File.join(@outdir, sam)
      cmd = "#{@bwa} mem"
      # options
      cmd << " -t #{@threads}"
      cmd << " #{index}"
      if right
        cmd << " #{left[:current]}"
        cmd << " #{right[:current]}"
      else
        cmd << " #{left[:current]}"
      end
      cmd << " > #{sam}"
      align = Cmd.new(cmd)
      align.run
      if align.status.success?
        left[:sam] = sam
        right[:sam] = nil if right
      else
        msg = "bwa failed\n#{align.stdout}\n#{align.stderr}"
        raise RuntimeError.new(msg)
      end
    end

    def get_stats
      str=""
      @stats.each do |key, value|
        str << "#{key}\t#{value}\n"
      end
      return str
    end

  end

end
