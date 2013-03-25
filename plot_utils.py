import matplotlib.pyplot as plt
import json
import argparse

mom_names = ["intensity", "mean", "sigma", "skewness"]


def plot_moments(prefix):
    data = json.load(open(prefix + ".json"))
    y = data["position"]
    # Intensities
    plt.subplot(411)
    plt.title("Velocity moments")
    plt.plot(y, data["intensity"], label='total')
    plt.plot(y, data["I1"], label='1')
    plt.plot(y, data["I2"], label='2')
    plt.ylim(0.0, None)
    plt.ylabel("Intensity")
    plt.grid()
    plt.legend(fontsize="x-small")
    plt.setp(plt.gca().get_xticklabels(), visible=False)

    # Velocities
    plt.subplot(412)
    plt.plot(y, data["mean"], label='mean')
    plt.plot(y, data["v1"], label='1')
    plt.plot(y, data["v2"], label='2')
    plt.ylabel("Velocity")
    plt.grid()
    plt.ylim(-5.0, 40.0)
    plt.legend(fontsize="x-small")
    plt.setp(plt.gca().get_xticklabels(), visible=False)

    # Sigma
    plt.subplot(413)
    plt.plot(y, data["sigma"])
    plt.ylabel("Sigma")
    plt.grid()
    plt.ylim(0.0, 16.0)
    plt.setp(plt.gca().get_xticklabels(), visible=False)

    # Skewness
    plt.subplot(414)
    plt.plot(y, data["skewness"])
    plt.ylabel("Skew")
    plt.grid()

    plt.xlabel("Position along slit, arcsec")
    plt.gcf().set_size_inches((4, 12))
    plt.savefig(prefix + "-plot.pdf", bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Exercise the plotting routines"""
    )
    parser.add_argument(
        "stampname", type=str,
        help="""Prefix of stamp file (e.g., p84-O_I_5577)"""
    )
    parser.add_argument(
        "--suffix", type=str, default="nc",
        help="""Suffix of stamp image"""
    )
    cmd_args = parser.parse_args()
    plot_moments("Stamps/{}-stamp-{}-moments".format(
        cmd_args.stampname, cmd_args.suffix))
