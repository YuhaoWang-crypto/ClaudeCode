# Modal from a sandboxed session — setup & troubleshooting

## Connecting through the agent proxy
Outbound HTTPS from the sandbox goes through a CONNECT proxy (`HTTPS_PROXY`). Modal
supports proxies natively but needs helper packages:

```bash
pip install modal 'python-socks[asyncio]' aiohttp-socks
```

| Symptom | Cause | Fix |
|---|---|---|
| `Could not connect to the Modal server` | grpclib can't use the proxy | `pip install 'python-socks[asyncio]'` |
| `ImportError: ... aiohttp-socks ...` when a call returns | large return value downloaded via aiohttp blob | `pip install aiohttp-socks` |
| Both above at once | — | `pip install 'modal[api-proxy-support]'` covers both |

`api.modal.com` is reachable through the proxy (curl returns 200); you do **not**
need any /etc/hosts hack or a manual TCP relay — Modal handles the proxy once the
socks helpers are installed.

## Auth
Credentials come from env: `MODAL_TOKEN_ID`, `MODAL_TOKEN_SECRET`. Verify with
`modal profile current` and `modal app list`.

## Running long jobs from the session
- `modal run` streams progress; pipe to a file and poll it, or run in the
  background. Do NOT wrap in `tail -N` mid-pipe — it buffers all output until EOF.
- Modal progress uses `\r`; reveal it with `tr '\r' '\n'`.
- ⚠️ **Ephemeral apps keep running server-side even if the local client dies.** If you
  kill/background a `modal run` and relaunch, old apps may crash-loop and warm
  containers can serve **stale mounted code**. Symptom: a fixed file still errors on
  the old line. Fix: `modal app stop <ap-id> -y` for every stale app (list with
  `modal app list`), and/or bump the `modal.App("name")` to force fresh containers.
- ⚠️ Beware `pkill -f "<pattern>"` matching the `modal run --stage <pattern>` command
  line itself — it will kill your own client.

## Cost hygiene
Stop orphaned apps when done: `modal app list` then `modal app stop <ap-id> -y`.
GPU functions bill while a container is alive; a disconnected client can leave one
running.
