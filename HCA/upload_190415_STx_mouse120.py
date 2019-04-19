import os
from typing import Mapping, Tuple, Union
import numpy as np
from starfish.types import Axes, Coordinates, Features, Number
from starfish import Codebook
from starfish.experiment.builder import FetchedTile, TileFetcher
from slicedimage import ImageFormat
from skimage.io import imread
from starfish.experiment.builder import write_experiment_json

# imaging
channels = ["DAPI", "Cy3", "Cy5", "AF750", "AF488"]
FILENAME = ["HCA_09_ms_120_c1_stitched",
            "HCA_09_ms_120_c2_stitched",
            "HCA_09_ms_120_c3_stitched",
            "HCA_09_ms_120_c4_stitched",
            "HCA_09_ms_120_c5_stitched"]

# codebook in csv (column 1: code, column2-x: color code in each round)
# AVOID comma in gene names (column 1)
codebook_csv = r"E:\HCA_09b_mouse_120genes\codebook_starfish_mouse120.csv"
DO_decorators = ["AF750", "AF488", "Cy3", "Cy5"]
CHORDER = [channels.index(i) for i in DO_decorators]

# tile position and metadata
# TODO: read metadata using bio-format
pixelscale = 0.1625
tilepos_xy_csv = r"E:\HCA_09b_mouse_120genes\Preprocessing\ReslicedTiles\tilepos.csv"

# not sure why this is so important and can't be read from image itself...
tilesz = 2000
SHAPE = tilesz, tilesz

# file location
input_dir = r"E:\HCA_09b_mouse_120genes\Preprocessing\ReslicedTiles"
output_dir = r"E:\HCA_09b_mouse_120genes\upload_starfish_190415\ISS_m_brain_03"


# hakuna matata
if not os.path.isdir(output_dir):
    try:
        os.mkdir(output_dir)
    except:
        os.makedirs(output_dir)


def add_codebook(experiment_json_doc):
    experiment_json_doc['codebook'] = "codebook.json"
    return experiment_json_doc


def make_codebook_json(codebook_csv):
    """ convert color code matrix in csv to json format"""
    codebook_array = []
    with open(codebook_csv, "r") as f:
        for line in f:
            line = line.rstrip('\n').split(',')
            codewords = []
            for r, colorcode in enumerate(line[1:]):
                codewords.append({Axes.ROUND.value: r,
                                  Axes.CH.value: int(colorcode)-1,
                                  Features.CODE_VALUE:1})
            codebook_array.append({Features.CODEWORD:codewords, Features.TARGET: line[0]})
    codebook = Codebook.from_code_array(codebook_array)
    codebook_json_filename = "codebook.json"
    codebook.to_json(os.path.join(output_dir, codebook_json_filename))


def get_tilepos(tilepos_xy_csv):
    tilexy = np.ndarray((0,2))
    with open(tilepos_xy_csv, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split(',')
            tilexy = np.vstack([tilexy, [np.double(line[0]), np.double(line[1])]])
    return tilexy


# TODO: skip reading and writing of images
class ISSTile2D(FetchedTile):
    def __init__(self, file_path, fov):
        self.file_path = file_path
        self.fov = fov

    @property
    def shape(self) -> Tuple[int, ...]:
        return SHAPE

    @property
    def coordinates(self) -> Mapping[Union[str, Coordinates], Union[Number, Tuple[Number, Number]]]:
        return {
            Coordinates.X: (tilexy[self.fov, 0]*pixelscale, (tilexy[self.fov, 0] + tilesz)*pixelscale),
            Coordinates.Y: (tilexy[self.fov, 1]*pixelscale, (tilexy[self.fov, 1] + tilesz)*pixelscale),
            Coordinates.Z: (0.0, 0.0),
        }

    def tile_data(self) -> np.ndarray:
        return imread(self.file_path)


class ISS2DPrimaryTileFetcher(TileFetcher):
    def __init__(self, path):
        self.path = path

    def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
        return ISSTile2D(os.path.join(self.path, "{}-{}/tile{}.tif".format(FILENAME[r], CHORDER[ch]+1, fov+1)), fov)


class ISS2DAuxTileFetcher(TileFetcher):
    def __init__(self, path, filename_prefix):
        self.path = path
        self.prefix = filename_prefix

    def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
        return ISSTile2D(os.path.join(self.path, self.prefix + "{}.tif".format(fov+1)), fov)


tilexy = get_tilepos(tilepos_xy_csv)


write_experiment_json(
    path=output_dir, fov_count=704, tile_format=ImageFormat.TIFF,
    primary_image_dimensions={
        Axes.ROUND: 5,
        Axes.CH: 4,
        Axes.ZPLANE: 1,
    },
    aux_name_to_dimensions={
        'nuclei': {
            Axes.ROUND: 1,
            Axes.CH: 1,
            Axes.ZPLANE: 1,
        },
        # 'dots': {
        #     Axes.ROUND: 1,
        #     Axes.CH: 1,
        #     Axes.ZPLANE: 1,
        # },
    },
    primary_tile_fetcher=ISS2DPrimaryTileFetcher(input_dir),
    aux_tile_fetcher={
        'nuclei': ISS2DAuxTileFetcher(input_dir, FILENAME[0] + "-1/tile"),
        # 'dots': ISS2DAuxTileFetcher(input_dir, FILENAME[0] + "-5/tile"),
    },
    postprocess_func=add_codebook,
    default_shape=SHAPE
)

make_codebook_json(codebook_csv)

