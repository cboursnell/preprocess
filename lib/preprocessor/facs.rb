require 'fileutils'

module Preprocessor

  class Facs

    def initialize(outdir, threads, memory, filter,
                   k, false_positive, threshold)
      @outdir = outdir
      @threads = threads
      @memory = memory
      @filter = filter # file or url of reference sequences
      @k = k
      @false_positive = false_positive
      @threshold = threshold
      @gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      @facs = File.join(@gem_dir, 'bin', 'facs')
      if !File.exist?(@facs)
        install
      end
      if File.exist?(File.expand_path(@filter))
        @filter = File.expand_path(@filter)
        # file
      else
        # try url
        file = "#{@outdir}/#{File.basename(@filter)}"
        cmd = "wget #{@filter} -O #{file}"
        dl = Cmd.new(cmd)
        dl.run
        if !dl.status.success?
          raise RuntimeError.new("Couldn't download #{@filter}")
        end
        @filter = file
      end
      # check if a file called @filter exists
      #   if it does use it to build the index
      #   if it doesn't try to download it
      #   if a filter wasn't specified use the default rrna fasta file
    end

    def install
      # download zip filehttps://github.com/SciLifeLab/facs/archive/master.zip
      dl = "#{@gem_dir}/bin/facs.zip"
      cmd = "wget https://github.com/SciLifeLab/facs/archive/master.zip"
      cmd << " -O #{dl}"
      if !File.exist?(dl)
        wget = Cmd.new(cmd)
        wget.run
      end
      # unpack zip
      cmd = "unzip -qq #{dl} -d #{@gem_dir}/bin/"
      unzip = Cmd.new(cmd)
      unzip.run
      # run make to compile
      Dir.chdir "#{@gem_dir}/bin/facs-master" do
        make = Cmd.new("make")
        make.run
      end
      # copy out facs binary
      cmd = "cp #{@gem_dir}/bin/facs-master/facs/facs #{@gem_dir}/bin/."
      copy = Cmd.new(cmd)
      copy.run
      # delete facs source code
      File.delete(dl) if File.exist?(dl)
      cmd = "rm -rf #{@gem_dir}/bin/facs-master"
      remove_dir = Cmd.new(cmd)
      remove_dir.run
    end

    def run left, right=nil
      # check for bloom filter index
      bloom = "#{@outdir}/#{File.basename(@filter)}.bloom"
      if !File.exist?(bloom)
        # build index if it doesn't exist
        cmd = "#{@gem_dir}/bin/facs build"
        cmd << " -r #{@filter}"
        cmd << " -o #{bloom}"
        cmd << " -e #{@false_positive}"
        build = Cmd.new(cmd)
        build.run
      end

      # run facs remove using reads and filter
      if left
        cmd = "#{@gem_dir}/bin/facs remove"
        cmd << " -r #{bloom}"
        cmd << " -q #{left[:current]}"
        cmd << " -t #{@threshold}"
        cmd << " -o #{@outdir}/"
        output_file = "#{@outdir}/"
        output_file << "#{File.basename(left[:current]).split(".").first}"
        output_file << "_#{File.basename(@filter)}"
        output_file << "_clean.fastq"
        remove = Cmd.new(cmd)
        remove.run
        left[:current] = output_file
      end
      if right
        cmd = "#{@gem_dir}/bin/facs remove"
        cmd << " -r #{bloom}"
        cmd << " -q #{right[:current]}"
        cmd << " -t #{@threshold}"
        cmd << " -k #{@k}" if @k
        cmd << " -o #{@outdir}/"
        output_file = "#{@outdir}/"
        output_file << "#{File.basename(right[:current]).split(".").first}"
        output_file << "_#{File.basename(@filter)}"
        output_file << "_clean.fastq"
        remove = Cmd.new(cmd)
        remove.run
        right[:current] = output_file
      end
    end

  end

end