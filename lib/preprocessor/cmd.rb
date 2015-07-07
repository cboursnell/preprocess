require 'open3'

module Preprocessor

  class Cmd

    attr_accessor :cmd, :stdout, :stderr, :status

    def initialize cmd
      @cmd = cmd
    end

    def run file=nil
      unless file.nil?
        @stdout, @stderr, @status = Open3.capture3 "echo #{file} exists"
        return true if File.exist?(file)
      end
      @stdout, @stderr, @status = Open3.capture3 @cmd
      return false
    end

    def to_s
      @cmd
    end

  end

end