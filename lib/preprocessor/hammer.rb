require 'yaml'

module Preprocessor

  class Hammer

    def initialize(outdir, threads, memory)
      @outdir = outdir
      @threads = threads
      @memory = memory

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      spades = File.join(gem_dir, 'bin', 'SPAdes-3.1.0-Linux',
                         'bin', 'spades.py')
      if !File.exist?(spades)
        wget_cmd = "wget http://spades.bioinf.spbau.ru/release3.1.0/"
        wget_cmd << "SPAdes-3.1.0-Linux.tar.gz -O "
        wget_cmd << "#{File.join(gem_dir,"bin","spades.tar.gz")}"

        tar_cmd = "tar xzf #{File.join(gem_dir,"bin","spades.tar.gz")}"
        tar_cmd << " -C #{File.join(gem_dir,"bin")}"
        download = Cmd.new(wget_cmd)
        download.run
        if download.status.success?
          # puts download.cmd
          # puts download.stdout

          unpack = Cmd.new(tar_cmd)
          unpack.run
          if unpack.status.success?
            # puts unpack.cmd
            # puts unpack.stdout
          end
        else
          puts download.stdout
          puts download.stderr
        end
      end
    end

    def run left, right=nil
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      spades = File.join(gem_dir, 'bin', 'SPAdes-3.1.0-Linux',
                         'bin', 'spades.py')
      dir = File.join(@outdir,
                      "hammer-#{left[:name]}-#{left[:type]}-#{left[:rep]}")
      cmd = "#{spades} -1 #{left[:current]} -2 #{right[:current]}"
      cmd << " --pe1-s #{left[:unpaired]}" if left[:unpaired]
      cmd << " --pe2-s #{right[:unpaired]}" if right[:unpaired]
      cmd << " --only-error-correction"
      cmd << " --disable-gzip-output"
      cmd << " -t #{@threads}"
      cmd << " -m #{@memory}"
      cmd << " -o #{dir}"
      # puts cmd
      hammer_cmd = Cmd.new(cmd)
      hammer_cmd.run
      if !hammer_cmd.status.success?
        raise RuntimeError.new("BayesHammer failed\n#{hammer_cmd.stderr}")
      end
      yaml = YAML.load_file("#{dir}/corrected/corrected.yaml")
      cat_cmd = "cat "
      left[:current] = yaml[0]["left reads"][0]
      right[:current] = yaml[0]["right reads"][0]

      if yaml[0] and yaml[0]["single reads"]
        yaml[0]["single reads"].each do |single|
          cat_cmd << " #{single} "
        end
        out = "#{dir}/left_unpaired.fq"
        cat_cmd << " >  #{out}"
        cat = Cmd.new(cat_cmd)
        cat.run
        left[:unpaired] = "#{out}"
      end
      if yaml[1] and yaml[1]["single reads"]
        yaml[1]["single reads"].each do |single|
          cat_cmd << " #{single} "
        end
        out = "#{dir}/right_unpaired.fq"
        cat_cmd << " >  #{out}"
        cat = Cmd.new(cat_cmd)
        cat.run
        right[:unpaired] = "#{out}"
      end

    end

  end

end
