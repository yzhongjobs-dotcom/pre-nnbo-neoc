"""Read existing source PDFs; keep page anchors and provenance for the audit."""
from pathlib import Path
import hashlib
import json
import fitz

SOURCE_ROOT = Path(r"D:/桌面/科研Agent Harness设计/impedance-control-e2e/source-materials")
OUTPUT = Path(__file__).resolve().parent
SOURCES = {
    "beber": "user-batch-mpvic/A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography.pdf",
    "pollayil": "user-batch-mpvic/Choosing_Stiffness_and_Damping_for_Optimal_Impedance_Planning.pdf",
    "wijayarathne": "user-batch-mpvic/Real-Time_Deformable-Contact-Aware_Model_Predictive_Control_for_Force-Modulated_Manipulation.pdf",
    "haninger": "user-batch-mpvic/Model predictive impedance control with Gaussian processes for human and environment interaction.pdf",
    "xue": "user-batch-mpvic/Model predictive variable impedance control towards safe robotic interaction in unknown disturbance-rich environments.pdf",
    "risk_sensitive": "user-batch-mpvic/Impedance Optimization for Uncertain Contact Interactions Through Risk Sensitive Optimal Control.pdf",
    "fu": "user-batch-mpvic/Human-Inspired_Active_Compliant_and_Passive_Shared_Control_Framework_for_Robotic_Contact-Rich_Tasks_in_Medical_Applications.pdf",
    "piedra_volumetric": "user-batch-topvenue-calibration/Visuotactile and Explicitly Force-Controlled Robotic Ultrasound for Abdominal Volumetric Reconstruction.pdf",
}

def main():
    out = OUTPUT / "原文定位"
    out.mkdir(exist_ok=True)
    manifest = []
    for key, relative in SOURCES.items():
        path = SOURCE_ROOT / relative
        with fitz.open(path) as doc:
            pages = [page.get_text(sort=True) for page in doc]
            (out / f"{key}.txt").write_text("\n\n".join(
                f"===== PDF PAGE {i + 1} =====\n{body}" for i, body in enumerate(pages)
            ), encoding="utf-8")
            manifest.append({"key": key, "source": str(path), "pages": len(pages),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "extraction": str(out / f"{key}.txt")})
            print(f"{key}: {len(pages)} pages, {sum(map(len, pages))} characters")
    (OUTPUT / "local_source_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
