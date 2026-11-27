import matplotlib.pyplot as plt
from wfn.propagation import PeResults


fname = "./results/wape/t0/100/side_0000.h5"
res = PeResults.from_file(fname)
ax = res.plot_side(freq=1000,vmin=-80, vmax=-20, add_colorbar=False)
plt.plot(res.r, res.s, 'k', lw=0.8)


fname = "./results/wape/t0/100/side_0179.h5"
res = PeResults.from_file(fname)
res.ds["r"] = -res.ds["r"]
res.plot_side(freq=1000,vmin=-80, vmax=-20,
                   cbar_kwargs={"label":"$\Delta L$ (dB)", "shrink": 0.3})
plt.plot(res.r, res.s, "k", lw=0.8)
plt.xlim(-500,500)
plt.gca().set_aspect("equal")
plt.tight_layout()

plt.show()
