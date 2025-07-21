from enum import Enum


class SpectralSystemEnum(str, Enum):
    C2SWAN = "C2_swan"
    N2CB = "N2CB"
    N2PlusBX = "N2PlusBX"
    NHAX = "NHAX"
    NOBX = "NOBX"
    OHAX = "OHAX"

    @property
    def display_name(self) -> str:
        """Human-readable display name for this spectral system."""
        display_names = {
            "C2_swan": "C₂ (d-a)",
            "N2CB": "N₂ (C-B)",
            "N2PlusBX": "N₂⁺ (B-X)",
            "NHAX": "NH (A-X)",
            "NOBX": "NO (B-X)",
            "OHAX": "OH (A-X)",
        }
        return display_names.get(self.value, self.value)

    @property
    def filename(self) -> str:
        """Database filename for this spectral system."""
        return f"{self.value}.db"
