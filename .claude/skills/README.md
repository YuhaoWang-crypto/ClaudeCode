# Project skills

Skills in this directory are auto-discovered by Claude Code for this repository.

## Vendored zine / poster design skills

Four third-party image-design skills were vendored here from public GitHub repos.
Each skill directory name matches its `SKILL.md` frontmatter `name`, which is what
Claude Code requires for discovery.

| Skill directory | Upstream repo | License |
| --- | --- | --- |
| `scene-distillation-zine-v1-3` | [Zeejay0/gathered-scenes-zine-skill](https://github.com/Zeejay0/gathered-scenes-zine-skill) | Gathered Scenes Zine Personal Non-Commercial License v1.0 |
| `scenes-gathered-zine-v1-3` | [Zeejay0/gathered-scenes-zine-skill](https://github.com/Zeejay0/gathered-scenes-zine-skill) | Gathered Scenes Zine Personal Non-Commercial License v1.0 |
| `photo-abstract-editorial` | [ZzzLc0405/photo-abstract-editorial](https://github.com/ZzzLc0405/photo-abstract-editorial) | Personal / non-commercial use only |
| `gc-minimal-zine-poster-v0-3` | [LiamGvchi/gc-minimal-zine-poster](https://github.com/LiamGvchi/gc-minimal-zine-poster) | MIT |

**Two of these are non-commercial-only.** The upstream `LICENSE` file is kept inside
each skill directory; read it before using the skill for anything commercial.

### What was vendored

Everything a skill needs at runtime: `SKILL.md`, `references/`, `agents/`, and
`evals/` where present, plus `assets/examples/` for `photo-abstract-editorial`
because its `SKILL.md` links to it.

Repo-level material that no `SKILL.md` reads was left upstream: the README files,
the example galleries in `gathered-scenes-zine-skill/examples/` and
`gc-minimal-zine-poster/examples/`, the brand images, and the donation QR codes.
Follow the upstream links above to see those.

### Updating

Re-clone the upstream repo and copy the same paths over the skill directory.
