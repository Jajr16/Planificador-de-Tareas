"""
URL configuration for django_bio project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
# from django.contrib import admin
from django.urls import path
from initial_setup import views

# Define all urls through the file views
urlpatterns = [
    # path('admin/', admin.site.urls),
    path('', views.config, name="home"),
    path('Sudoku', views.sudoku_start, name="Sudoku"),
    path('UpSudo', views.sudoku_Update, name="sudoku_Update"),
    path('Sim', views.start_simulation, name="Sudoku_Simulation"),
    path('get_csrf_token/', views.get_csrf_token, name='get_csrf_token'),
]
