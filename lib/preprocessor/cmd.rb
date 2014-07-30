require 'open3'

module Preprocessor

  class Cmd

    attr_accessor :cmd, :stdout, :stderr, :status

    def initialize cmd
      @cmd = cmd
    end

    def run
      @stdout, @stderr, @status = Open3.capture3 @cmd
    end

  end

end