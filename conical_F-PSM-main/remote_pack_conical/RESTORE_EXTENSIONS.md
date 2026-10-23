# RESTORE_EXTENSIONS.md - read this first, then restore the file names

*Written 2026-10-23 on the benchmark PC. Every file in this archive except this one has had `.txt` appended to
its name so that e-mail accepts the archive. Nothing else was changed: the SHA-256 in the table is of each
file's contents, which renaming does not touch. British English.*

## Agent instructions

1. Work inside the extracted `remote_pack_conical` folder; this file sits at its root.
2. Run ONE of the two scripts below from that folder. Each reads the table at the end of this file, renames
   every file back to its original name, checks its SHA-256, prints one line per file and a final count.
   Running it twice is harmless: a file already restored is only re-checked.
3. The run must end with `restored 63 OK, 0 BAD`. If any file is `BAD` or `MISSING`, stop and report it;
   do not continue with a damaged pack.
4. Then open `CLAUDE.md` (now restored) and follow it.

### Python (3.8 or later)

Save as `restore.py` beside this file and run `python restore.py`:

```python
import os, hashlib, re
ok = bad = 0
for line in open("RESTORE_EXTENSIONS.md", encoding="utf-8"):
    m = re.match(r"^\| `([^`]+)` \| `([^`]+)` \| `[^`]*` \| `([0-9a-f]{64})` \|", line)
    if not m:
        continue
    src, dst, want = m.groups()
    if os.path.isfile(src):
        os.replace(src, dst)
    if not os.path.isfile(dst):
        print("MISSING", dst); bad += 1; continue
    got = hashlib.sha256(open(dst, "rb").read()).hexdigest()
    if got == want:
        ok += 1; print("OK ", dst)
    else:
        bad += 1; print("BAD", dst)
print(f"restored {ok} OK, {bad} BAD")
```

### PowerShell

Save as `restore.ps1` beside this file and run `powershell -ExecutionPolicy Bypass -File restore.ps1`:

```powershell
$ok = 0; $bad = 0
Get-Content RESTORE_EXTENSIONS.md -Encoding UTF8 | ForEach-Object {
  if ($_ -match '^\| `([^`]+)` \| `([^`]+)` \| `[^`]*` \| `([0-9a-f]{64})` \|') {
    $src = $Matches[1]; $dst = $Matches[2]; $want = $Matches[3]
    if (Test-Path -LiteralPath $src) { Move-Item -LiteralPath $src -Destination $dst -Force }
    if (-not (Test-Path -LiteralPath $dst)) { Write-Output "MISSING $dst"; $bad++; return }
    $got = (Get-FileHash -LiteralPath $dst -Algorithm SHA256).Hash.ToLower()
    if ($got -eq $want) { $ok++; Write-Output "OK  $dst" } else { $bad++; Write-Output "BAD $dst" }
  }
}
Write-Output "restored $ok OK, $bad BAD"
```

## Table: 63 files

Paths are relative to the pack root, forward slashes. `original extension` is informational; the
scripts use `original name`.

| file in this archive | original name | original extension | sha256 |
|---|---|---|---|
| `CLAUDE.md.txt` | `CLAUDE.md` | `.md` | `ca38ab2f623fda5071fda3be578079ad1f7aac61924b4c42eae1ff07e9ed31e2` |
| `EXECUTION_CONICAL_LUMERICAL.md.txt` | `EXECUTION_CONICAL_LUMERICAL.md` | `.md` | `d3a56c8cde20b304afda0f28ca5ae1d589bc7e9f2734546d0ca14716f690260e` |
| `README_PACK.md.txt` | `README_PACK.md` | `.md` | `43ef9adeb69c1e4cf03461f5c7294f94e0f567e51c1b00eeee95b4af5c4c3dbd` |
| `RETURN_CONICAL.md.txt` | `RETURN_CONICAL.md` | `.md` | `5e82375cfa7b0cb03187af8c4ec4f0657a36f284cb010a377b629918b886188e` |
| `code/basis.py.txt` | `code/basis.py` | `.py` | `843b57c7a4640a4ef8be3bb8bf2e5cbc7131fbddcb306b7a59142d30375a4742` |
| `code/build_pack.py.txt` | `code/build_pack.py` | `.py` | `0e3937da015debf7be82d0c67c40b974b9660263660d4a4c94aa807b6ac957ff` |
| `code/cases.py.txt` | `code/cases.py` | `.py` | `d0038831b29b5958f5c6267df3379ab4af01e36dd1deebb505f4d72b91e13445` |
| `code/cases_conical.py.txt` | `code/cases_conical.py` | `.py` | `a7fae2f2d19a87581d4bc1dd1e2453fb39e96f3332214a491f19d3e2caaff096` |
| `code/check_env.py.txt` | `code/check_env.py` | `.py` | `da8412f89efdeb1fd1a696b711095e17239d28af946bb037532b181c38751eb6` |
| `code/compare.py.txt` | `code/compare.py` | `.py` | `0e3253deb92c30c805946ae2ff8a91a6054e68a575b7926714ec851c1fd6bc8c` |
| `code/compare_conical.py.txt` | `code/compare_conical.py` | `.py` | `d94e363926dc6575187348cb4717ad7392eb4895c1ca2cd9e22db5b08dac7550` |
| `code/dryrun_conical.py.txt` | `code/dryrun_conical.py` | `.py` | `429302e7c22eb45ba632216fc2be3c7ff37a132d90adf6b62b24c5592842bfcd` |
| `code/figures_conical.py.txt` | `code/figures_conical.py` | `.py` | `e685ce1a1cbc61b24313e349f2a2821d486edc2c71f07b2666fccaf1c4b9db72` |
| `code/fpsm.py.txt` | `code/fpsm.py` | `.py` | `459b9cc807069853d972c9b392f51c17a75019e1f048e86dcc5024ac8f5363c8` |
| `code/fpsm_conical.py.txt` | `code/fpsm_conical.py` | `.py` | `27fae60d6db9c10380e70d1080a8b810cb89e3766c08771504ea41906321a55b` |
| `code/kogelnik.py.txt` | `code/kogelnik.py` | `.py` | `506fd9ff3af155f705da7c78d0f6b21abf7fcb1a58d3ff8cf2ed2b459fe7b00c` |
| `code/kogelnik_te.py.txt` | `code/kogelnik_te.py` | `.py` | `0a566c637c685082beda305c49c0d17860f6b45444e7ecd0a7e869bd59f4493b` |
| `code/lumerical_rcwa.py.txt` | `code/lumerical_rcwa.py` | `.py` | `0d9ec221bd5dc36bce05c1a796ad4a12b67b01175ea7941ee5db2dc2172ba945` |
| `code/lumerical_rcwa_conical.py.txt` | `code/lumerical_rcwa_conical.py` | `.py` | `f88c58ca5f0d388964fd6124cd2de0608980f55a0e0e4e5ac77a2ec518ad97bb` |
| `code/map2d.py.txt` | `code/map2d.py` | `.py` | `24f8855b7bd30f342a15f89d4848f75acbbdceeae4a21732034c4e95144fbb02` |
| `code/mock_lumerical.py.txt` | `code/mock_lumerical.py` | `.py` | `fd6389569d2ce46374e5eeacf5e0e1a381f28ae8fd14d9d219e629946c82d7c8` |
| `code/rcwa.py.txt` | `code/rcwa.py` | `.py` | `f689bf6da241fe31816ca94132f18d848cb360083fe3a72881e9433d4b88eba3` |
| `code/rcwa_conical.py.txt` | `code/rcwa_conical.py` | `.py` | `3a9b100421b64e40759f236540053579fce749e9a49ca71154f44b2cde0c038c` |
| `code/rcwa_te.py.txt` | `code/rcwa_te.py` | `.py` | `4dfe4c4da9e4a3da085949e5b413fad37dbef9732d9d9818fb56a124f62b8f0c` |
| `code/reference_values_tm.py.txt` | `code/reference_values_tm.py` | `.py` | `08b91f2e3f6030597bf03a64c9fd62ac4aac1593e49f6c6998a41ce41e862603` |
| `code/stress_scan.py.txt` | `code/stress_scan.py` | `.py` | `ee9f5a49906aa298593bce780ab226fc161400c38978ca502e4d49e6cf852b58` |
| `code/stress_summary.py.txt` | `code/stress_summary.py` | `.py` | `b47f6b77b6c1113a06b318e3a9cd7eb56193dd34fd2c96ec1b68a7ba6e6a527a` |
| `code/summary_chart.py.txt` | `code/summary_chart.py` | `.py` | `1818d80ede73b0fded6a6c3feae5e85286f44dbf55ee1b5138afbe3f2c35e141` |
| `code/sweep_david.py.txt` | `code/sweep_david.py` | `.py` | `a75c88692c8dbad0df04c9b0a22417ecfe8416933d705ab0d4c3bbd620b0632c` |
| `code/test_conical.py.txt` | `code/test_conical.py` | `.py` | `90ba8d7aefcbdfb6986c1ba03cbd9bf4076c40e3293fe100d37535621ec0be5d` |
| `code/test_engines.py.txt` | `code/test_engines.py` | `.py` | `d29babf3871de006d95a45f84bda242d656585827d34eaa94256c8593867324c` |
| `code/test_fpsm_vectors.py.txt` | `code/test_fpsm_vectors.py` | `.py` | `fc430ffe1489653519f9b4f3e162a00ec77c9ed1b0af336513ccd02ca9b47d5e` |
| `code/test_kogelnik_pol.py.txt` | `code/test_kogelnik_pol.py` | `.py` | `673e64bae047041be47e6c3688170be65bc1c975bb0a1623f6b5cf1c76116b02` |
| `code/test_tm_referee.py.txt` | `code/test_tm_referee.py` | `.py` | `4e8d85e9b286189de74080faf2462516b96de0f60624e2b75a0c1ed57e24386d` |
| `code/zemax_kogelnik_nsc.py.txt` | `code/zemax_kogelnik_nsc.py` | `.py` | `7b9dbe3acca26e98f2c641963ee0da84d80c7d7ab5c9cc3fedac9fa66da13eb8` |
| `code/zemax_rcwa_nsc.py.txt` | `code/zemax_rcwa_nsc.py` | `.py` | `be27eb74c80201b2b0eef1ea542b1f9088a4b8a6db0900327a0de9c45c222e69` |
| `code/zemax_rcwa_probe2.py.txt` | `code/zemax_rcwa_probe2.py` | `.py` | `b1ae6cfcca1c7b2a2e274f95511cabd736c09a94edc00c8d08a10928c9bd8c3b` |
| `code/zos_connect.py.txt` | `code/zos_connect.py` | `.py` | `1d4e54015ee799cf5b77c962a54d329ff1eb841d70929f0a7467b87a5862b858` |
| `evidence/conical_gates.txt.txt` | `evidence/conical_gates.txt` | `.txt` | `7afbee973c59c6631771a1be1f1e2ae98b0085a2d1ec15d75791f5fe92ddeb16` |
| `figures/reference_conical_bench.png.txt` | `figures/reference_conical_bench.png` | `.png` | `8630ef0938750c7f3433da24f74e0fe3473f2ddcf1e3bb28c304dc7e46685968` |
| `figures/sweep_errors_vs_mu.png.txt` | `figures/sweep_errors_vs_mu.png` | `.png` | `707e8846d4845bdb48f51c859525f2447bb5f17d14205260135102858b8691e8` |
| `lumerical/ansys_example/license.txt.txt` | `lumerical/ansys_example/license.txt` | `.txt` | `3539737e9548f41f5474612614246f7acbdee264c71805f37491f191a3aeeed4` |
| `lumerical/ansys_example/metadata.json.txt` | `lumerical/ansys_example/metadata.json` | `.json` | `5db9eb615cca74eac4b93039bd4e19ac57c67c6d59858753aa52e60c77ec352b` |
| `lumerical/ansys_example/vhg_R_6um_model.fsp.txt` | `lumerical/ansys_example/vhg_R_6um_model.fsp` | `.fsp` | `e593b7368ea326840700e08644da0201ad5075ad8e96d036c9abba4dc6acc7e5` |
| `reference/Conical_F_PSM_Mixed_Report.pdf.txt` | `reference/Conical_F_PSM_Mixed_Report.pdf` | `.pdf` | `74ebb297491eca855d0e1bf008c07c1c1f86f0901b73f0a11771602d0e4895d8` |
| `reference/Conical_F_PSM_RCWA_Report.pdf.txt` | `reference/Conical_F_PSM_RCWA_Report.pdf` | `.pdf` | `d3aba759d3403e46e6d7fc7a65ce88571b1cf6167c820395022dbf437cbb5456` |
| `reference/EXECUTION_LUMERICAL_original.md.txt` | `reference/EXECUTION_LUMERICAL_original.md` | `.md` | `9b6cfb2efb9bebcf624a15aed0398bd50909ea758196a6545aeadf2dfe95f289` |
| `reference/ansys_article_36965778575379_layer_repetition.txt.txt` | `reference/ansys_article_36965778575379_layer_repetition.txt` | `.txt` | `b99099fa928c91559357a26ceb7b5c23a730e9f6d17cc8f42b7100d1142fa639` |
| `reference/ansys_article_37118679173267_VHG_lumerical.txt.txt` | `reference/ansys_article_37118679173267_VHG_lumerical.txt` | `.txt` | `09dddfb691a533ec705a8d8a369e5ef7f55d1e38a5ec7aae33d9da2f82fe98f9` |
| `reference/ansys_article_42661816933139_zemax_rcwa_dll.txt.txt` | `reference/ansys_article_42661816933139_zemax_rcwa_dll.txt` | `.txt` | `455d0176608f2a0c31147243310e5e14fcb736ddc46ab8d4b84404cfa4ee4dc3` |
| `reference/ansys_rcwa_verdict_opticstudio.md.txt` | `reference/ansys_rcwa_verdict_opticstudio.md` | `.md` | `cd5e361d83be455faf7f54ce9745b7c5990d3490721a6d24ea729dac04262eb9` |
| `reference/execution_original.md.txt` | `reference/execution_original.md` | `.md` | `bd540381042abde512315615737f6471b21479ec93cf3c4f1fe2399eed18b464` |
| `reference/lumerical_converge_lumerical_R_2026-09-18.txt.txt` | `reference/lumerical_converge_lumerical_R_2026-09-18.txt` | `.txt` | `a41dd0f38e08c22e37cc0948d06bbd219c8e02703d8a68eeb2fe6a75a4f37ae7` |
| `reference/lumerical_probe_2026-09-18.txt.txt` | `reference/lumerical_probe_2026-09-18.txt` | `.txt` | `bdb8acb356cda0bb5a88dbf06d367e02e9012cdfb6305f64b1185a211c6cf3a0` |
| `reference/lumerical_verdict.md.txt` | `reference/lumerical_verdict.md` | `.md` | `1c868c87ec7ea6b31833e716b6660ab8532599087fd5e3ba111d4b9a94e2a360` |
| `reference/summary_2026-09-30.csv.txt` | `reference/summary_2026-09-30.csv` | `.csv` | `bb1e0770dd8b73b7ecd344b6729ce0427ecda3d12d993a9e83c1cb34f6502c04` |
| `reference/tm_pack_kogelnik_gates.txt.txt` | `reference/tm_pack_kogelnik_gates.txt` | `.txt` | `c4614edfa05dc9910c910890b6f6105f2078c5eb26a405efe056d8da7d8fb7cd` |
| `reference/tm_pack_tm_referee_gates.txt.txt` | `reference/tm_pack_tm_referee_gates.txt` | `.txt` | `38c108dbc40e4a1fbb321e3330a22487911d8f789916c6281bb8ceaf0647723e` |
| `reference/tm_pack_tm_reference_values.md.txt` | `reference/tm_pack_tm_reference_values.md` | `.md` | `eec1003ef91ed7602df55cb93a040cb31e8176f004e48a197dcf39955fc7e32f` |
| `reference/zemax_rcwa_probe2_2026-09-18.txt.txt` | `reference/zemax_rcwa_probe2_2026-09-18.txt` | `.txt` | `a935333d54f9f9bacd60158d866af867a840550fb6d9378cbefdb407dceb1230` |
| `results/mg_fig5_thickness_scan.csv.txt` | `results/mg_fig5_thickness_scan.csv` | `.csv` | `9a2c3affd9dc62214c3c46ee4c3495527e96802ef599f8a5fae431f0316c6ce1` |
| `results/reference_conical.csv.txt` | `results/reference_conical.csv` | `.csv` | `d5ddb29f1097eb73278003fb33c58e0092329d311ebfa35ac3cba8cae93fdcb1` |
| `results/sweep_504.csv.txt` | `results/sweep_504.csv` | `.csv` | `8ecdee2f705f5d12b40884c1eadc6b2f0eeaeb151579c3631c1dd9f58d162a0d` |
