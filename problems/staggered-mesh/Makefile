all:
	g++ main.cpp -O3
gif:
	convert frame/*.png anime.gif
movie:
	ffmpeg -r 10 -i frame/%04d.png -vcodec mjpeg -sameq output.avi
