---
license: apache-2.0
task_categories:
- robotics
tags:
- LeRobot
- ur
- robotics
configs:
- config_name: default
  data_files: data/*/*.parquet
---

This dataset was created using [LeRobot](https://github.com/huggingface/lerobot).

## Dataset Description



- **Homepage:** [More Information Needed]
- **Paper:** [More Information Needed]
- **License:** apache-2.0

## Dataset Structure

[meta/info.json](meta/info.json):
```json
{
    "codebase_version": "v2.1",
    "robot_type": "ur",
    "total_episodes": 90,
    "total_frames": 15527,
    "total_tasks": 5,
    "total_videos": 0,
    "total_chunks": 1,
    "chunks_size": 1000,
    "fps": 30,
    "splits": {
        "train": "0:90"
    },
    "data_path": "data/chunk-{episode_chunk:03d}/episode_{episode_index:06d}.parquet",
    "video_path": "videos/chunk-{episode_chunk:03d}/{video_key}/episode_{episode_index:06d}.mp4",
    "features": {
        "base_rgb": {
            "dtype": "image",
            "shape": [
                224,
                224,
                3
            ],
            "names": [
                "height",
                "width",
                "channel"
            ]
        },
        "wrist_rgb": {
            "dtype": "image",
            "shape": [
                224,
                224,
                3
            ],
            "names": [
                "height",
                "width",
                "channel"
            ]
        },
        "ultrasound_rgb": {
            "dtype": "image",
            "shape": [
                480,
                640,
                3
            ],
            "names": [
                "height",
                "width",
                "channel"
            ]
        },
        "state": {
            "dtype": "float32",
            "shape": [
                12
            ],
            "names": [
                "x",
                "y",
                "z",
                "roll",
                "pitch",
                "yaw",
                "fx",
                "fy",
                "fz",
                "tx",
                "ty",
                "tz"
            ]
        },
        "actions": {
            "dtype": "float32",
            "shape": [
                6
            ],
            "names": [
                "dx",
                "dy",
                "dz",
                "droll",
                "dpitch",
                "dyaw"
            ]
        },
        "prompt": {
            "dtype": "string",
            "shape": [
                1
            ],
            "names": [
                "prompt"
            ]
        },
        "us_gt": {
            "dtype": "int32",
            "shape": [
                1
            ],
            "names": [
                "us_gt"
            ]
        },
        "stage": {
            "dtype": "int32",
            "shape": [
                1
            ],
            "names": [
                "stage"
            ]
        },
        "timestamp": {
            "dtype": "float32",
            "shape": [
                1
            ],
            "names": null
        },
        "frame_index": {
            "dtype": "int64",
            "shape": [
                1
            ],
            "names": null
        },
        "episode_index": {
            "dtype": "int64",
            "shape": [
                1
            ],
            "names": null
        },
        "index": {
            "dtype": "int64",
            "shape": [
                1
            ],
            "names": null
        },
        "task_index": {
            "dtype": "int64",
            "shape": [
                1
            ],
            "names": null
        }
    }
}
```


## Citation

**BibTeX:**

```bibtex
[More Information Needed]
```