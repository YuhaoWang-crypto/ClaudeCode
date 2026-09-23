---
title: Laya API
emoji: 🎯
colorFrom: indigo
colorTo: blue
sdk: docker
app_port: 7860
startup_duration_timeout: 1h
pinned: false
license: apache-2.0
---

# Laya API

Serves [convaiinnovations/laya](https://huggingface.co/convaiinnovations/laya) over HTTP with `laya-serve`.

- `GET /health`
- `POST /v1/systemone` (requires `Authorization: Bearer <LAYA_API_KEY>` when that secret is set)

```bash
curl https://<user>-laya-api.hf.space/v1/systemone \
  -H "Authorization: Bearer $LAYA_API_KEY" \
  -H 'content-type: application/json' \
  -d '{"state":{"body":"billed twice, refund please"},
       "questions":{"dept":{"type":"choice","instructions":"which team?",
                    "criteria":{"billing":"refunds","tech":"bugs"}}}}'
```
