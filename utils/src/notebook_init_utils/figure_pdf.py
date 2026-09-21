import pathlib

from matplotlib.backends.backend_pdf import PdfPages


class FigurePDFs:
    """Append figures to one multi-page PDF per output path.

    Each PDF is written to `<name>.pdf.partial` and only renamed to its final
    name by `close()`, so an interrupted run never leaves a truncated,
    unopenable PDF behind (and `path.exists()` means the PDF is complete).
    Raster elements are embedded at `dpi` (default 600).
    """

    def __init__(self, dpi: int = 600):
        self.dpi = dpi
        self._pdfs: dict[pathlib.Path, PdfPages] = {}

    @staticmethod
    def _partial(pdf_path: pathlib.Path) -> pathlib.Path:
        return pdf_path.with_name(pdf_path.name + ".partial")

    def savefig(self, fig, pdf_path, dpi=None, **kwargs) -> None:
        pdf_path = pathlib.Path(pdf_path)
        if pdf_path not in self._pdfs:
            pdf_path.parent.mkdir(parents=True, exist_ok=True)
            self._pdfs[pdf_path] = PdfPages(self._partial(pdf_path))
        self._pdfs[pdf_path].savefig(fig, dpi=dpi or self.dpi, **kwargs)

    def close(self, pdf_path=None) -> None:
        """Finalise one PDF (or all of them when `pdf_path` is None)."""
        paths = list(self._pdfs) if pdf_path is None else [pathlib.Path(pdf_path)]
        for path in paths:
            pdf = self._pdfs.pop(path, None)
            if pdf is not None:
                pdf.close()
                self._partial(path).replace(path)
