#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <vmd_tachyon_trajectory_files_directory_path>"
  exit 1
fi

DIR=$1

if [ ! -d "$DIR" ]; then
  echo "Error: $DIR is not a directory"
  exit 1
fi

cd "$DIR" || exit

shopt -s extglob
if ! ls untitled.+([0-9]).ppm 1> /dev/null 2>&1; then
  echo "Error: No files matching pattern 'untitled.%05d.ppm' found in directory"
  exit 1
fi

# Create video from VMD movie maker (Renderer=Tachyon, Movie Settings=Trajectory) output via ffmpeg
ffmpeg -r 60 -i untitled.%05d.ppm -vcodec libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -crf 0 -pix_fmt yuv420p my_video.mp4

# Slow MP4 down via ffmpeg (4.0 == 4.0x times slower)
ffmpeg -i my_video.mp4 -vf "setpts=4.0*PTS" my_slowed_video.mp4

rm -f untitled.+([0-9]).ppm
rm -f untitled.+([0-9]).dat
rm -f untitled.par
