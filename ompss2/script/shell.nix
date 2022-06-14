let
  pkgs = import (builtins.fetchTarball
    "https://pm.bsc.es/gitlab/rarias/bscpkgs/-/archive/master/bscpkgs-master.tar.gz");

  rWrapper = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [ tidyverse rjson jsonlite egg viridis ];
  };

  # Recursively set MPI
  bsc' = pkgs.bsc.extend (self: super: {
    mpi = self.impi;
    #mpi = self.openmpi;
  });

#  clangOmpss2UnwrappedFixed = bsc'.clangOmpss2UnwrappedGit.overrideAttrs (old: rec {
#    src = builtins.fetchGit {
#      url = "ssh://git@bscpm03.bsc.es/llvm-ompss/llvm-mono.git";
#      ref = "master";
#
#      # Broken
#      #rev = "dc297872575e16afcd526120118a365bab150efc";
#
#      # Was working okeish
#      #rev = "ecc7282a0f7f8494366e42dbc710ffda388ccec9";
#
#      # 2022-04-14: Testing HEAD to see if I can avoid a crash
#      rev = "d4a6748b53036787166cd6957b5f1dd16f8379a5";
#    };
#    version = src.shortRev;
#  });
#
#  clangOmpss2Fixed = bsc'.clangOmpss2Git.override {
#    clangOmpss2Unwrapped = clangOmpss2UnwrappedFixed;
#  };

  extrae4 = bsc'.extrae.overrideAttrs (old: rec {
    version = "3.7.1";
    src = pkgs.fetchFromGitHub {
      owner = "bsc-performance-tools";
      repo = "extrae";
      rev = "${version}";
      sha256 = "sha256-aoGM8yRE3KBDHEZOzPmIQzIzCWcencWTYWV00jRPKsw=";
    };
  });

  hwloc-1-11-6 = bsc'.callPackage ~/bscpkgs/bsc/hwloc/1.11.6/default.nix {};
  slurm-16-05-8-1 = bsc'.callPackage ~/bscpkgs/bsc/slurm/16.05.8.1/default.nix {
    hwloc = hwloc-1-11-6;
  };

  oldslurm = slurm-16-05-8-1;

  gaspi = /nix/store/j01fzm5i5w6f0zhdxbwfkw7f173rv061-GPI-2-f5eb152;
  tagaspi = /nix/store/x01f2gml9k7pmhkibrxj05hdhpas2lf3-tagaspi-5aabb18;

  clangOmpss2Fixed =
    /nix/store/hj2ig1dhahjw1mcc3l385ky26hvalfz8-clang-ompss2-wrapper-d4a6748;
in
  pkgs.mkShell {
    name = "minimd";
    NIX_HARDENING_ENABLE = "";
    buildInputs = with bsc'; [ pkgs.python3 babeltrace2 nanos6
    extrae4 mpi icc
    mcxx
    clangOmpss2Fixed pkgs.cmake
    llvmPackages.lldb
    pkgs.gdb
    rWrapper
    oldslurm
    tampi
    gaspi
    tagaspi
  ] ++ (with pkgs; [ vim ]);
    shellHook = ''
      echo "NOTE: using mpi=${bsc'.mpi}"
    '';
  }
