module Preprocessor

  class Rcorrector

    def initialize(outdir, threads, memory)
      @outdir = outdir
      @threads = threads
      @memory = memory

      @rcorrector = Which.which('run_rcorrector.pl').first
      if @rcorrector.nil?
        abort "Couldn't find Rcorrector in PATH"
      end
      @kmer = 23
    end

    def run left, right
      leftfiles = []
      left.each do |data|
        leftfiles << File.expand_path(data[:current])
      end
      rightfiles = []
      right.each do |data|
        rightfiles << File.expand_path(data[:current])
      end
      cmd = "#{@rcorrector} "
      cmd << " -1 #{leftfiles.join(",")}"
      cmd << " -2 #{rightfiles.join(",")}"
      cmd << " -k #{@kmer}"
      cmd << " -t #{@threads}"
      cmd << " -od #{@outdir}"

      # TODO check if output files already exist
      correct = Cmd.new cmd
      missing = false
      (leftfiles+rightfiles).each do |file|
        file = File.basename(file, File.extname(file))
        file = File.join(@outdir, file + ".cor.fq")
        unless File.exist?(file)
          missing = true
        end
      end
      if missing
        correct.run
        unless correct.status.success?
          puts "Something went wrong running Rcorrector"
          puts correct.stdout
          puts correct.stderr
          abort "oh well"
        end
      else
        puts "Corrected files already exist"
      end
      # TODO set output files as current
      (left+right).each do |data|
        file = data[:current]
        file = File.basename(file, File.extname(file))
        file = File.join(@outdir, file + ".cor.fq")
        data[:current] = file
      end
    end

	end

end


# Usage: perl ./run_rcorrector.pl [OPTIONS]
# OPTIONS:
# Required parameters:
#         -s seq_files: comma separated files for single-end data sets
#         -1 seq_files_left: comma separated files for the first mate in the paried-end data sets
#         -2 seq_files_right: comma separated files for the second mate in the paired-end data sets
#         -i seq_files_interleaved: comma sperated files for interleaved paired-end data sets
# Other parameters:
#         -k kmer_length (default: 23)
#         -od output_file_directory (default: ./)
#         -t number of threads to use (default: 1)
#         -maxcorK INT: the maximum number of correction within k-bp window (default: 4)
#         -ek expected_number_of_kmers: does not affect the correctness of program but affect the memory usage (default: 100000000)      -stdout: output the corrected reads to stdout (default: not used)
#         -stage INT: start from which stage (default: 0)
#                 0-start from begining(storing kmers in bloom filter);
#                 1-start from count kmers showed up in bloom filter;
#                 2-start from dumping kmer counts into a jf_dump file;
#                 3-start from error correction.
