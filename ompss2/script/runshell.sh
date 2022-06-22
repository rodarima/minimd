#!/bin/sh

export NIX_SSL_CERT_FILE=/etc/ssl/ca-bundle.pem
dir=$(dirname $0)

nix-setup nix-shell --pure -v --tarball-ttl 9999999 "$dir/shell.nix"
