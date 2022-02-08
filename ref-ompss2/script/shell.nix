let
  pkgs = import (builtins.fetchTarball
    "https://pm.bsc.es/gitlab/rarias/bscpkgs/-/archive/nosv/bscpkgs-master.tar.gz");

  rWrapper = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [ tidyverse rjson jsonlite egg ];
  };
in
  pkgs.mkShell {
    name = "minimd";
    buildInputs = with pkgs.bsc; [ nanos6 extrae impi icc mcxx tampi rWrapper ];
    shellHook = ''
      export LANG=C
    '';
  }
