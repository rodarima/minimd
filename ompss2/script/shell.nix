let
  pkgs = import (builtins.fetchTarball
    "https://pm.bsc.es/gitlab/rarias/bscpkgs/-/archive/update-ompss2/bscpkgs-master.tar.gz");

  rWrapper = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [ tidyverse rjson jsonlite egg viridis ];
  };

  mpi = pkgs.bsc.impi;
  #mpi = pkgs.bsc.openmpi;

  clangOmpss2UnwrappedFixed = pkgs.bsc.clangOmpss2UnwrappedGit.overrideAttrs (old: rec {
    src = builtins.fetchGit {
      url = "ssh://git@bscpm03.bsc.es/llvm-ompss/llvm-mono.git";
      ref = "master";
      #rev = "dc297872575e16afcd526120118a365bab150efc";
      rev = "ecc7282a0f7f8494366e42dbc710ffda388ccec9";
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
    buildInputs = with pkgs.bsc; [ pkgs.python3 babeltrace2 nanos6
    extrae mpi icc
    mcxx
    clangOmpss2Fixed pkgs.cmake
    rWrapper (tampi.override {mpi=mpi;}) ];
    shellHook = ''
      export LANG=C
    '';
  }
