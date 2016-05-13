require 'fixwhich'
require 'bindeps'

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

      @java = Which.which('java').first
      raise RuntimeError.new("java not installed") if !@java
      version = Cmd.new("java -version")
      version.run
      version_number = 0
      if version.stderr =~ /([0-9._\-]+)/
        version_number = $1
        version_number = version_number.split(".")[0..1].join(".").to_f
      end
      if version_number and version_number < 1.7
        msg = "bbnorm requires java version 1.7 or higher\n"
        msg << "You have #{version_number}"
        raise RuntimeError.new(msg)
      end

      @gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      @bbnorm = File.join(@gem_dir, 'bin', 'bbmap', 'bbnorm.sh')
      if !File.exist?(@bbnorm)
        install
      end
    end

    def install
      url = "http://downloads.sourceforge.net/project/bbmap/BBMap_35.43.tar.gz"
      dl = "#{@gem_dir}/bin/bbmap.tar.gz"
      cmd = "wget #{url} -O #{dl}"
      wget = Cmd.new(cmd)
      wget.run
      if wget.status.success?
        cmd = "tar xzf #{dl} --directory #{@gem_dir}/bin"
        untar = Cmd.new(cmd)
        untar.run
        if !untar.status.success?
          raise RuntimeError.new("untar failed")
        end
        File.delete(dl)
      else
        raise RuntimeError.new("download failed")
      end
    end

    def run left, right=nil
      if left and right
        cmd = "#{@bbnorm} "
        cmd << "in=#{left[:current]} "
        cmd << "in2=#{right[:current]} "
        leftout = "#{left[:name]}-"
        leftout << "#{left[:type]}-" if left[:name]!=left[:type]
        leftout << "#{left[:rep]}_#{left[:pair]}.bbnorm.fq"
        rightout = "#{right[:name]}-"
        rightout << "#{right[:type]}-" if right[:name]!=right[:type]
        rightout << "#{right[:rep]}_#{right[:pair]}.bbnorm.fq"
        leftout = File.join(@outdir, leftout)
        rightout = File.join(@outdir, rightout)
        cmd << "out=#{leftout} "
        cmd << "out2=#{rightout} "
        cmd << "k=#{@k} "
        cmd << "hashes=#{@tables} "
        cmd << "bits=#{@bits} "
        cmd << "lowthresh=#{@lowthresh} "
        cmd << "mindepth=#{@mindepth} "
        cmd << "minkmers=#{@minkmers} "
        cmd << "target=#{@coverage} "
        cmd << "threads=#{@threads}"
        norm = Cmd.new(cmd)
        unless File.exist?(leftout)
          norm.run
          if !norm.status.success?
            msg = "something went wrong with bbnorm\n"
            msg << "#{norm.stdout}\n#{norm.stderr}"
            raise RuntimeError.new(msg)
          end
        end
        left[:prenorm] = left[:current]
        right[:prenorm] = right[:current]
        left[:current] = leftout
        right[:current] = rightout
      end
    end

  end

end
