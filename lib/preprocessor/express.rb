require 'bindeps'
require 'csv'
require 'fixwhich'

module Preprocessor

  class Express


    def initialize(outdir, threads, reference)
      @outdir = outdir
      @threads = threads
      @reference = File.expand_path(reference)

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'express.yaml')
      Bindeps.require gem_deps
      @express = Which.which('express').first
    end

    def run left, right=nil
      sam = left[:sam]
      express_out = File.join(@outdir, "express_#{File.basename(left[:file])}")
      cmd = "#{@express}"
      cmd << " --output-dir #{express_out}"
      cmd << " --no-update-check "
      cmd << " --additional-batch 2 "
      cmd << " --output-align-prob "
      cmd << " #{@reference}"
      cmd << " #{sam}"
      express = Cmd.new(cmd)
      express.run
      if !express.status.success?
        msg = "Something went wrong with eXpress\n#{express.stdout}\n"
        msg << "#{express.stderr}"
        raise RuntimeError.new(msg)
      end
      left[:express] = File.join(express_out, "results.xprs")
    end

    def get_stats
      str=""
      return str
    end

  end

end
