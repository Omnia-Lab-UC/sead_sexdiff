{
  description = "sead-sexdiff project development environment flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
  };

  outputs =
    {
      self,
      nixpkgs,
      systems,
    }:
    let
      forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in
    {
      devShells = forEachSystem (
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          default = pkgs.mkShell {
            packages = with pkgs; [
              python312
              uv
              R
              rPackages.SingleCellExperiment
              rPackages.MAST
              rPackages.anndataR
              rPackages.edgeR
              rPackages.limma
              rPackages.DESeq2
            ];

            LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath [
              pkgs.stdenv.cc.cc.lib
              pkgs.zlib
            ];

            shellHook = ''
              export UV_PYTHON="${pkgs.python312}/bin/python"
              export UV_PYTHON_DOWNLOADS=never

              uv sync --locked
            '';
          };
        }
      );
    };
}
