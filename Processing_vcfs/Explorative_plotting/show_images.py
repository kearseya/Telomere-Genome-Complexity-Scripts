#!/usr/bin/python3
import click
import matplotlib.pyplot as plt
import pandas as pd
import cv2
import os

@click.command()
@click.argument("indir", default="/path/to/project/variants/processing_vcfs/circos_plots/icgc")
@click.option("-p", "--prediction", default="/path/to/project/variants/icgc/allpredot.csv")
@click.option("-r", "--rows", default=4)
@click.option("-c", "--cols", default=5)
def show_split(indir, prediction, rows, cols):
	indexes = {"ssidx": [], "lsidx": [], "osidx": []}
	images_d = {"s": [], "l": [], "o": []}
	png_list = os.listdir(indir)
	sample_list = [i.split("_")[0] for i in png_list]
	if prediction != None:
		split_df = pd.read_csv(prediction)
		short_samples = split_df[split_df["short"] == True]["sample"].tolist()
		long_samples = split_df[split_df["short"] == False]["sample"].tolist()
		#print("s", short_samples, "l", long_samples)
		for i, s in enumerate(sample_list):
			if s in short_samples:
				indexes["ssidx"].append(i)
			elif s in long_samples:
				indexes["lsidx"].append(i)
			else:
				indexes["osidx"].append(i)
		for t in ["ssidx", "lsidx", "osidx"]:
			for i in indexes[t]:
				images_d[t[0]].append(cv2.imread(os.path.join(indir, png_list[i])))
		print(indexes)
	else:
		images_d["o"] = [cv2.imread(os.path.join(indir, i)) for i in png_list]

	for z in images_d:
		images = images_d[z]
		print("len = ", len(images))
		for i in range(0, len(images), rows*cols):
			print("i = ", i)
			fig = plt.figure(figsize=(8,8))
			fig.tight_layout()
			if z == "s":
				fig.suptitle("Short")
			elif z == "l":
				fig.suptitle("Long")
			for j in range(0, rows*cols):
				print("j = ", j)
				fig.add_subplot(rows, cols, j+1)
				if i+j+1 <= len(images):
					plt.imshow(images[i+j])
					#plt.tight_layout()
					plt.title(sample_list[i+j][:3])
				plt.axis("off")
			plt.show()


if __name__ == "__main__":
	show_split()

