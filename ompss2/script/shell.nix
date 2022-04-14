let
  pkgs = import (builtins.fetchTarball
    "https://pm.bsc.es/gitlab/rarias/bscpkgs/-/archive/master/bscpkgs-master.tar.gz");

  rWrapper = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [ tidyverse rjson jsonlite egg viridis ];
  };

  #mpi = pkgs.bsc.impi;
  mpi = pkgs.bsc.openmpi;

  clangOmpss2UnwrappedFixed = pkgs.bsc.clangOmpss2UnwrappedGit.overrideAttrs (old: rec {
    src = builtins.fetchGit {
      url = "ssh://git@bscpm03.bsc.es/llvm-ompss/llvm-mono.git";
      ref = "master";

      # Broken
      #rev = "dc297872575e16afcd526120118a365bab150efc";

      # Was working okeish
      #rev = "ecc7282a0f7f8494366e42dbc710ffda388ccec9";

      # 2022-04-14: Testing HEAD to see if I can avoid a crash
      rev = "d4a6748b53036787166cd6957b5f1dd16f8379a5";
    };
    version = src.shortRev;
  });

  extrae4 = pkgs.bsc.extrae.overrideAttrs (old: rec {
    version = "3.7.1";
    src = pkgs.fetchFromGitHub {
      owner = "bsc-performance-tools";
      repo = "extrae";
      rev = "${version}";
      sha256 = "sha256-aoGM8yRE3KBDHEZOzPmIQzIzCWcencWTYWV00jRPKsw=";
    };
  });

  clangOmpss2Fixed = pkgs.bsc.clangOmpss2Git.override {
    clangOmpss2Unwrapped = clangOmpss2UnwrappedFixed;
  };
in
  pkgs.mkShell {
    name = "minimd";
    NIX_HARDENING_ENABLE = "";
    buildInputs = with pkgs.bsc; [ pkgs.python3 babeltrace2 nanos6
    extrae4 mpi icc
    mcxx
    clangOmpss2Fixed pkgs.cmake
    tagaspi
    gaspi gpi-2
    pkgs.gdb
    rWrapper (tampi.override {mpi=mpi;}) ];
    shellHook = ''
      echo "NOTE: using mpi=${mpi}"
    '';
  }
