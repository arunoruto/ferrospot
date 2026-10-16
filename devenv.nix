{
  pkgs,
  config,
  ...
}:
{
  env = {
    UV_PYTHON = toString config.languages.python.package.interpreter;
  };

  # pdfimages and pdftoppm, for ferrospot.pdf.
  packages = [ pkgs.poppler-utils ];

  scripts = {
    dashboard.exec = ''uv run ferrospot-dashboard "$@"'';
  };

  enterShell = ''
    git --version
    if [ ! -L "$DEVENV_ROOT/.venv" ]; then
        ln -s "$DEVENV_STATE/venv/" "$DEVENV_ROOT/.venv"
    fi
  '';

  languages.python = {
    enable = true;

    uv = {
      enable = true;
      sync = {
        enable = true;
        groups = [
          "test"
          "lint"
        ];
        extras = [
          "dashboard"
          "simulate"
        ];
      };
    };

    libraries = with pkgs; [
      zlib
      libGL
    ];
  };

  git-hooks.hooks = {
    ruff-format = {
      enable = true;
      description = "Ruff formatter";
      package = config.git-hooks.tools.ruff;
      entry = "ruff format";
      types = [ "python" ];
      args = [ "--check" ];
    };

    ruff-check = {
      enable = true;
      description = "Ruff linter";
      package = config.git-hooks.tools.ruff;
      entry = "ruff check";
      types = [ "python" ];
      after = [ "ruff-format" ];
    };
  };
}
