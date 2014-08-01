require 'bindeps'
require 'csv'

module Preprocessor

  class Bowtie2

    require 'which'
    include Which

    def initialize(outdir, threads, reference, expression)
      @outdir = outdir
      @threads = threads
      @reference = File.expand_path(reference)
      @expression = expression

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'bowtie2.yaml')
      Bindeps.require gem_deps
      @build = which('bowtie2-build').first
      @bowtie = which('bowtie2').first
      @stats={:reads=>0,:paired_reads=>0,:unaligned=>0,:unique=>0,:multi=>0}
    end

    def run left, right=nil
      # construct index if it doesn't exist
      index = File.join(@outdir, File.basename(@reference))
      if !File.exist?("#{index}.1.bt2")
        cmd = "#{@build}"
        cmd << " #{@reference}"
        cmd << " #{index}"
        build = Cmd.new(cmd)
        build.run
        if !build.status.success?
          msg = "Something went wrong with bowtie2-build\n"
          msg << "#{build.stdout}"
          msg << "#{build.stderr}"
          raise RuntimeError.new(msg)
        end
      end
      l = File.basename(left[:current])
      r = File.basename(right[:current])
      sam = "#{l.split(".")[0..-2].join(".")}-"
      sam << "#{r.split(".")[0..-2].join(".")}-"
      sam << "#{File.basename(index)}.sam"
      sam = File.join(@outdir, sam)
      cmd = "#{@bowtie}"
      cmd << " -x #{index}"
      cmd << " -1 #{left[:current]}"
      cmd << " -2 #{right[:current]}"
      cmd << " -p #{@threads}"
      cmd << " --very-sensitive "
      cmd << " --no-unal" # temp TODO remove
      if @expression # if this sam file is going to be used for expression later
        cmd << " -a -X 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4"
        cmd << " --no-discordant --no-mixed "
      end
      cmd << " -S #{sam}"
      align = Cmd.new(cmd)
      align.run
      hash = {}
      if align.status.success?
        out = align.stderr.split("\n")
        out.each do |row|
          if row =~ /([0-9]+)\sreads;\ of\ these:/
            @stats[:reads]+=$1.to_i
          end
          if row =~ /([0-9]+) \([0-9]+\.[0-9]+\%\) were\ paired;\ of\ these:/
            @stats[:paired_reads]+=$1.to_i
          end
          if row =~ /([0-9]+) \([0-9]+\.[0-9]+\%\)\ aligned\ concordantly\ 0\ times/
            @stats[:unaligned]+=$1.to_i
          end
          if row =~ /([0-9]+) \([0-9]+\.[0-9]+\%\)\ aligned\ concordantly\ exactly\ 1\ time/
            @stats[:unique]+=$1.to_i
          end
          if row =~ /([0-9]+) \([0-9]+\.[0-9]+\%\)\ aligned\ concordantly\ >1\ times/
            @stats[:multi]+=$1.to_i
          end
        end
      else
        msg = "Bowtie2 failed\n#{align.stdout}\n#{align.stderr}"
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
