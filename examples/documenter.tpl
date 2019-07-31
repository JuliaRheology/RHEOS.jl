{% extends 'markdown.tpl'%}


{% block stream %}
```output
{{ output.text | trim }}
```
{% endblock stream %}


{% block data_html scoped %}
```@raw html
{{ output.data['text/html'] }}
```
{% endblock data_html %}


{% block data_text scoped %}
```output
{{ output.data['text/plain'] | trim }}
```
{% endblock data_text %}
