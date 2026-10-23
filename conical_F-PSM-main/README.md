# conical_F-PSM

Private transfer repository: task packs built on the benchmark PC for the Ansys work PC. Keep it
private — the packs carry Ansys knowledge-base article texts, an Ansys example model and
unpublished manuscripts.

| archive | built | contents | sha256 |
|---|---|---|---|
| `remote_pack_conical_2026-10-23_txt.zip` | 2026-10-23 | the conical F-PSM Lumerical leg (63 files), every file renamed `<name>.txt`; `RESTORE_EXTENSIONS.md` at the pack root restores the names | `859f081e5a2b3916bd5c6ca9c971cc726b63b6e85de9afeb5e511cc05e0f9efb` |

## On the work PC

1. Clone or download, check the archive's SHA-256 against the table, extract it into a folder of
   its own.
2. Open `remote_pack_conical/RESTORE_EXTENSIONS.md` and run one of its two restore scripts from
   that folder; it must end with `restored 63 OK, 0 BAD`.
3. Open the restored `CLAUDE.md` and follow it (about three hours of Lumerical solver time).
