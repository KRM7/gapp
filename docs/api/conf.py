project = "gapp"
copyright = "2023, Krisztián Rugási"

highlight_language = "c++"
primary_domain = "cpp"
maximum_signature_line_length = 80

extensions = [ "breathe" ]

breathe_projects = {"gapp":"doxygen-out"}
breathe_default_project = "gapp"

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/KRM7/genetic-algorithms",
    "use_repository_button": True,
    "use_fullscreen_button": False,
    "use_download_button": False,
    "show_prev_next": False,
}

html_static_path = ["_resources"]
html_css_files = ["custom.css"]