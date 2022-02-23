let
  pkgs = import (builtins.fetchTarball
    "https://pm.bsc.es/gitlab/rarias/bscpkgs/-/archive/update-ompss2/bscpkgs-master.tar.gz");

  rWrapper = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [ tidyverse rjson jsonlite egg ];
  };

  clangOmpss2UnwrappedFixed = pkgs.bsc.clangOmpss2UnwrappedGit.overrideAttrs (old: rec {
    src = builtins.fetchGit {
      url = "ssh://git@bscpm03.bsc.es/llvm-ompss/llvm-mono.git";
      ref = "master";
      #rev = "dc297872575e16afcd526120118a365bab150efc";
    };
    version = src.shortRev;
  });

  clangOmpss2Fixed = pkgs.bsc.clangOmpss2Git.override {
    clangOmpss2Unwrapped = clangOmpss2UnwrappedFixed;
  };
in
  pkgs.mkShell {
    name = "minimd";
    NIX_HARDENING_ENABLE = "";
    buildInputs = with pkgs.bsc; [ nanos6 extrae openmpi icc mcxx
    clangOmpss2Fixed pkgs.cmake
    rWrapper (tampi.override {mpi=openmpi;}) ];
    shellHook = ''
      export LANG=C
    '';
  }
