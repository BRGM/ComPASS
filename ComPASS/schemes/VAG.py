from collections import namedtuple

from .._kernel import get_kernel

VAGVolumeDistribution = namedtuple(
    "VAGVolumeDistribution", ["cell", "fracture"], defaults=[0.075, 0.15]
)


class VAGScheme:
    def __init__(self, omega_Darcy={}, omega_Fourier={}):
        self.omega_Darcy = VAGVolumeDistribution(**omega_Darcy)
        self.omega_Fourier = VAGVolumeDistribution(**omega_Fourier)

    def __str__(self):
        return f"VAG scheme with omega Darcy = {self.omega_Darcy}, omega Fourier = {self.omega_Fourier}"

    def compute_transmissivities(self):
        kernel = get_kernel()
        kernel.VAGFrac_TransDarcy()
        if kernel.has_energy_transfer_enabled():
            kernel.VAGFrac_TransFourier()

    def compute_volumes(self):
        """
        VAG volumes are computed locally on each proc.
        There is no need to synchronise volume of ghost control volumes.
        cf. https://gitlab.inria.fr/compass/v4/ComPASS/-/issues/301
        """
        kernel = get_kernel()
        kernel.VAGFrac_VolsDarcy(self.omega_Darcy.cell, self.omega_Darcy.fracture)
        if kernel.has_energy_transfer_enabled():
            kernel.VAGFrac_VolsFourier(
                self.omega_Fourier.cell, self.omega_Fourier.fracture
            )
