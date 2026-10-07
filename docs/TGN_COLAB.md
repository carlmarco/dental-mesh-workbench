# Scoring ToothGroupNetwork with this project's metrics (D81)

ToothGroupNetwork (Lim et al., winner of 3DTeethSeg'22, https://github.com/limhoyeon/ToothGroupNetwork) needs a
CUDA GPU with at least 11 GB, so it runs on Colab, not on the Mac. Its repository has no license file: use it to
evaluate only, and don't redistribute its code or checkpoints. Teeth3DS is CC BY-NC-ND 4.0: upload the scans
privately (your own Drive or Colab session), never to a public location.

## Which scans

The challenge checkpoints were trained on the 3DTeethSeg'22 training lists, which differ from the Teeth3DS
split this project uses. Only scans in `private-testing-set.txt` are held out for the network:

| Part | Scans | Of these, in the challenge private test set |
|---|---|---|
| 6 (fresh; never used for development) | 300 | 84 |
| 5 (test split, already evaluated 5 times) | 300 | 90 |

The primary comparison is the 84 part-6 scans, which neither method has seen.

## 1. Package the scans (Mac)

```bash
python3 tools/tgn/prepare_input.py data/data_part_6 data/tgn_input_part6.zip
```

Upload `data/tgn_input_part6.zip` to your Google Drive.

## 2. Run the network (Colab, GPU runtime: Runtime > Change runtime type > T4 or better)

```bash
!git clone -b challenge_branch https://github.com/limhoyeon/ToothGroupNetwork.git
%cd ToothGroupNetwork
!pip install wandb open3d multimethod termcolor trimesh easydict gdown
!pip install --ignore-installed PyYAML
!cd external_libs/pointops && python setup.py install
!gdown --folder https://drive.google.com/drive/folders/15oP0CZM_O_-Bir18VbSM8wRUEzoyLXby
```

Unzip `ckpts(challenge).zip` into `ckpts/`, mount Drive, unzip the scans, and run the final-phase model:

```bash
from google.colab import drive; drive.mount('/content/drive')
!unzip -q /content/drive/MyDrive/tgn_input_part6.zip -d /content/scans
!python inference_final.py --input_path /content/scans/input --save_path /content/pred
!cd /content && zip -qr tgn_pred_part6.zip pred
```

Download `tgn_pred_part6.zip`.

Known risks, untested here:
- The repo was tested on PyTorch 1.7.1 / CUDA 11.0. If `pointops` fails to compile against Colab's newer PyTorch,
  the usual fix is removing the deprecated `THC` includes in `external_libs/pointops/src`. Otherwise install an
  older PyTorch that matches Colab's CUDA.
- All 84 part-6 cases contain only one jaw (the other jaw of each patient is not in part 6). If
  `inference_final.py` expects both jaws, patch its loop to skip the missing one.
- The README asks for the Y axis pointing backwards. Teeth3DS scans come in the challenge's frame, so no change
  should be needed. If the predictions look like noise, check the axes first.

## 3. Score (Mac)

```bash
unzip -q tgn_pred_part6.zip -d data/tgn_pred_part6
./build-release/tools/margin_eval score data/data_part_6 data/tgn_pred_part6/pred "ToothGroupNetwork (final_v1)"
```

Scans without a prediction are skipped, so the comparison runs on exactly the scans the network was given.
The output adds paired per-scan statistics against this project's graph cut on those same scans.
