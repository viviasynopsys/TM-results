import glob, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
src = r"C:\Users\vivia\code_base\17-09 lumerical pack\conical_F-PSM-main\remote_pack_mixed\figures"
out = r"C:\Users\vivia\.copilot\session-state\c9e05458-8f88-43c9-acce-1719857e2a40\files"
files = sorted(glob.glob(os.path.join(src, "mixed_[A-D][0-9][0-9]_*.png")))
print(len(files))
per = 12
for s in range(0, len(files), per):
    chunk = files[s:s + per]
    fig, axs = plt.subplots(4, 3, figsize=(18, 17))
    for ax in axs.ravel():
        ax.axis("off")
    for ax, f in zip(axs.ravel(), chunk):
        ax.imshow(mpimg.imread(f))
        ax.set_title(os.path.basename(f)[6:12], fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"sheet_{s // per + 1}.png"), dpi=70)
    plt.close(fig)
