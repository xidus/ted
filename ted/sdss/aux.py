#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sun 29 Sep 2013
#   Initial build.
#

import os

def HTML_create_cutout_display(cutout_files, path=None, summary=None):
    """Creates HTML document that displays the cutout images."""
    if path is None:
        return

    HTML_head = '''\
<!doctype html </>
<html>
<head>
<title>Cutouts</title>
<style rel="stylesheet" type="text/css">
span {
    display: block;
    float: left;
    border: 1px #ddd solid;
    background: #f6f6f6;
    padding: 10px;
    margin: 10px 10px;
    text-align: center;
}
div.clear {
    clear: both;
}
span img {
    border: 5px #fff solid;
    margin: 0;
    padding: 0;
}
div {
    text-align: center;
}
div span img {
    width: auto;
}
pre {
    display:block;
}
</style>
</head>
<body>
'''
    HTML_foot = '''\
</body>
</html>
'''
    HTML_body = ''

    if summary is not None:
        HTML_body += '\n<pre>{}</pre>'.format(summary)

    fname_stats = 'stats.png'
    ifname_stats = os.path.join(path, fname_stats)
    if os.path.isfile(ifname_stats):
        HTML_body += '''

<div>
  <span><img src="{0}" alt="{1}" title="{1}" /></span>
</div>
<div class="clear"></div>'''.format(fname_stats, '')

    for fname in cutout_files:
        HTML_body += '\n\n<span>\n<img src="{0}" alt="{1}" title="{1}" /><br />\n<em>{2}</em>\n</span>'.format(
            fname, os.path.splitext(fname)[0], fname.split('T')[0]
        )

    HTML_content = HTML_head + HTML_body + HTML_foot
    ofname_html = os.path.join(path, 'index.html')
    with open(ofname_html, 'w+') as fsock:
        fsock.write(HTML_content)
    return os.path.isfile(ofname_html)

