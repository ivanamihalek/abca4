
from PIL import Image
import os

#########################################

def convert_to_jpg(fnm):
	# note e includes "."
	f, e = os.path.splitext(fnm)
	new_fnm = f + ".jpg"
	try:
		Image.open(fnm).convert('RGB').save(new_fnm)
	except IOError as err:
		print("cannot convert", fnm, err)
		return None
	os.remove(fnm)
	return


def convert_frames(frames_home, subdir):
	outdir = "{}/{}".format(frames_home, subdir)
	for path, subdirs, files in os.walk(outdir):
		print(files[:5])
		for fnm in files:
			if fnm[-3:]=="png":
				convert_to_jpg("{}/{}".format(path,fnm))
	return


def subdir_prep(frames_home, subdir):
	outdir = "{}/{}".format(frames_home, subdir)
	#if os.path.exists(outdir): shutil.rmtree(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

