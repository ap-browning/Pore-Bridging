# Plot options
plt.style.use("default")

matplotlib.rc("axes", linewidth=0.5)
matplotlib.rc("axes", axisbelow=true)
matplotlib.rc("axes", labelcolor="#545454")
matplotlib.rc("axes", edgecolor="#aaaaaa")
matplotlib.rc("axes", labelsize=10.0)
matplotlib.rc("axes", grid=true)
matplotlib.rc("axes.formatter", limits=[-3,4])
matplotlib.rc("xtick", color="#000000")
matplotlib.rc("ytick", color="#000000")
matplotlib.rc("xtick.major", width=0.5)
matplotlib.rc("ytick.major", width=0.5)
matplotlib.rc("xtick.major", size=2.0)
matplotlib.rc("ytick.major", size=2.0)
matplotlib.rc("xtick.major", pad=2.0)
matplotlib.rc("ytick.major", pad=2.0)
matplotlib.rc("grid",   linewidth=0.5)
matplotlib.rc("grid",   color="#fafafa")
matplotlib.rc("font",   size=8.0)
matplotlib.rc("font",   family="Helvetica")
matplotlib.rc("pdf",    fonttype=42)
matplotlib.rc("savefig",transparent=true)
matplotlib.rc("legend", fontsize=8)
matplotlib.rc("legend", fancybox=false)
matplotlib.rc("legend", edgecolor="#545454")
matplotlib.rc("legend", frameon=false)



# Colours and styles to use
PoreSizeCols = ["#0FA88C","#1F77B4","#D62728","#FF7F0E"] #BB74C2
PoreSizeColsLight = ["#00C7B5","#50BDF5","#D66C58","#FFB175"] #BB74C2
PoreSizeSymb = [".","s","d","*"]
PoreSizeMkSz = [10,4,5,5]
Snms      = Dict(S .=> ["Density","Coverage","Circularity","Edge Density"])
Slim      = [[0.0,0.005],[0.0,1.0],[0.0,1.0],[0.0,0.005]]
Stks      = Dict(S .=> [collect(range(Slim[i]...,length=[6,5,5,6][i])) for i = 1:4])
Slim      = Dict(S .=> [diff(x)[1]*[-0.05,0.05] + x for x ∈ Slim])


# For subplot numbering
alphabet = "abcdefghijklmnopqrstuvwxyz"
function NumberPlots!(axs)
    for i = 1:prod(size(axs))
        permutedims(axs)[i].set_title(join(["(",alphabet[i],")"]),loc="left")
    end
end
